# Adjacency and boundary conditions

This chapter describes how to define bonds between sites in a lattice and how to use them to define boundary conditions.

## General bonds

There are several ways to define bonds between sites in a lattice. The most general way is to define a function that takes two sites and returns if they are connected. 

Generally, all objects that define bonds are subtypes of the `AbstractBonds` type. Here is a short overview of types of bonds that are supported by this package:
- [`AdjacencyMatrix`](@ref) — a simple matrix that defines if two sites are connected.
- [`SiteDistance`](@ref) — connectivity based on the distance between sites.
- [`NearestNeighbor`](@ref) — connects nearest neighbors on a given lattice. You can specify the order of the neighbors (e. g. first, second, etc.) by passing an integer to the constructor: `NearestNeighbor(2)`.
- [`Translation`](@ref) — connects two sites that are shifted by a certain vector.
- [`BravaisTranslation`](@ref) — connects two sites in a Bravais lattice unit cell. Unlike the `Translation`, it is defined in terms of the lattice axes and the site indices.

Let's begin with the [`AdjacencyMatrix`](@ref). It is a simple wrapper around a boolean-valued matrix that defines if two sites are connected. Here is an example of how to use it:

```@example 1
using LatticeModels, Plots
l = SquareLattice(4, 4)

# Define a nearest neighbor adjacency matrix
adj = AdjacencyMatrix(l, NearestNeighbor())

# Let's disconnect the center sites from the rest vertically
for x in 2:3
    site1 = l[x=x, y=2]
    site2 = l[x=x, y=1]
    adj[site1, site2] = false

    site3 = l[x=x, y=3]
    site4 = l[x=x, y=4]
    adj[site3, site4] = false
end

plot(adj)                   # Plot what we've got
plot!(l, showbonds=false)   # And the lattice itself
```

As you can see, the adjacency matrix is a writable object, so you can change the bonds as you like.
In this example we deleted vertical bonds between the four center sites and the boundary of the lattice.

There are many ways to create an adjacency matrix — for example, you can use a function that takes two sites and returns if they are connected. This is useful when you need to define bonds in a more complex way:

```@example 1
adj2 = AdjacencyMatrix(l) do site1, site2
    dx = site1.x - site2.x
    dy = site1.y - site2.y
    return abs(dx * dy) == 1    # Diagonal bonds
end
plot(adj2)
plot!(l, showbonds=false)
```

## Translations

A "translation" is a special type of bonds that has a distinct direction (e. g. each pair of sites consists of a "source" and a "target" site). Also it must define exactly one destination for each source site (if any).

An example is the [`Translation`](@ref) — it allows defining bonds between sites that are shifted by a certain vector. 

```@example 2
using LatticeModels, Plots
l = GrapheneRibbon(6, 4)        # A convenient constructor for a honeycomb lattice
tr = Translation(l, [1, 2√3/3])
plot(tr)
plot!(l, ls=:dash, linecolor=:grey)
```

Another type of translation is the [`BravaisTranslation`](@ref). It translates sites on a Bravais lattice in terms of the lattice vectors:

```juila
tr2 = BravaisTranslation(1 => 2, [0, 1])
tr1 == tr2  # true
```

What does this notation mean? We take the first site in the unit cell; then we go to the unit cell shifted by `[0, 1]` and take the second site from there. Note that the first example described the same translation, but in terms of positions, not unit cell indices.

!!! note
    You can omit the pair of indices if you want to translate the unit cell regardless of the site index. 
    `Bravais[j1, j2, ...]` is shorthand for `BravaisTranslation([j1, j2, ...])`.    

You can use translations to shift a site, for example: `site2 = site1 + tr`, just like with regular vectors. Another 
use-case for translations is in defining boundary conditions, as we will see in the next section.

## [Boundary conditions](@id BoundaryConditions_chapter)

The most general form of boundary conditions supported by this package is this:

```math
\psi(r + R) = \psi(r) f(r)
```

where `R` is a translation vector and `f` is some function. This is a generalization of periodic boundary conditions, which are a special case of this form.

These conditions are in fact applied not to the lattice itself, but to the Hamiltonian. Think of it as a way to
replace the ``c(r)^\dagger c(r')`` hopping with ``c(r)^\dagger c(r' - R) f(r')``, if the ``r'`` site is not present in the lattice, but ``r' - R`` is.

In general, there are three types of boundary conditions:
- [`PeriodicBoundary`](@ref) — the most common type of boundary conditions. Just periodicity with no factor.
- [`TwistedBoundary`](@ref) — periodicity with a phase factor that does not depend on ``r``.
- [`FunctionBoundary`](@ref) —  general form of boundary conditions.

The constructor for the boundary accepts two arguments — the phase (or the function) and the translation. 

As an example, let us consider the example with the Hofsadter butterfly from the [Examples](@ref) page. The magnetic field `B` changes from zero to one flux quantum per plaquette, for each value of `B` we calculate the energy spectrum of an infinite lattice and plot it. 

Let's consider a periodic lattice instead. We want to apply magnetic field in the Landau gauge to it, and since 
the translation operators include the vector potential ``\mathcal{A}``, we have to tweak the boundary conditions a little:

```math
\psi(x + L_x, y) = \psi(x, y) e^{2\pi i B y L_x},
\psi(x, y + L_y) = \psi(x, y)
```

This can be done with the `FunctionBoundary`:

```@example 3
using LatticeModels, Plots
l = HoneycombLattice(10, 10)
Lx = 10
Ly = 5√3
n_plaquettes = 100
B_step = 1 / (Lx * Ly)          # Field step: one flux quantum through all plaquettes
B_max = n_plaquettes * B_step   # Until one flux quantum per plaquette

points_E = Float64[]
points_B = Float64[]
for B in 0:B_step:B_max
    f(site) = exp(2π * im * B * site.y * Lx)
    xboundary = FunctionBoundary(f, [Lx, 0])
    yboundary = PeriodicBoundary([0, Ly])
    lb = setboundaries(l, xboundary, yboundary)
    H = tightbinding_hamiltonian(lb, field=LandauGauge(B))
    dg = diagonalize(H)
    append!(points_E, dg.values)
    append!(points_B, fill(B, length(dg.values)))
end
scatter(points_B, points_E, xlabel="B", ylabel="E", leg=false, ms=1)
```

Note that we could have set the boundary conditions to the lattice in one line. In fact, we already did this in the
Examples section:

```julia
lb = setboundaries(l, [Lx, 0] => f, [0, Ly] => true)
```

Let's explain this notation:
- Each pair consists of a translation on the left and a "boundary specifier" on the right. A simple vector, like here, is interpreted as a [`Translation`](@ref) — you can use any translation type here, `Bravais[-5, 10]` instead of `[0, Ly]` would work as well.
- The second argument can be one of the following:
    - `true/false` — a periodic or open boundary condition.
    - A number `θ` — a twisted boundary condition with a phase factor `exp(im * θ)`.
    - A function `f` — a general boundary condition.
- Sometimes the lattice includes default translation axes, and you can use their aliases as translations. 
  For example, for a 10x10 square lattice `setboundaries(l, :axis1 => f, :axis2 => true)` is equivalent to `setboundaries(l, [10, 0] => f, [0, 10] => true)`.

Note that you can pass boundary conditions to the lattice constructor as well. We could not do this in the previous example, because the boundary conditions depended on the magnetic field `B`. However, if you have a fixed boundary condition, you can pass it to the lattice constructor like this:

```julia
lb = HoneycombLattice(10, 10, boundaries=([Lx, 0] => f, [0, Ly] => true))
```

This is equivalent to the `setboundaries` call from the previous example. 

!!! note
    In our example we imposed periodic boundary conditions in shape of a rectangle ``0 < x < L_x``, ``0 < y < L_y``.
    However, the lattice itself is more like a parallelogram, because of its unit cell shape. This is not a problem,
    because this shape is still periodic in terms of these translations along x and y.

    If the periodicity is violated, for example, if the both `r` and `r + R` sites are present in the lattice, an error will be thrown.

## Generic lattice

The `GenericLattice` is a tool you might want to use if you need to define a lattice with a more complex geometry. It is basically just a collection of arbitrary sites. As an example, let's define a lattice with a naive random geometry:

```@example 4
using LatticeModels, Plots
l = GenericLattice{2}()         # 2D lattice
for i in 1:100
    pt = rand(2) * 5
    md = minimum(site -> norm(pt - site.coords), l, init=1.0)
    md > 0.5 && push!(l, pt)    # Add a site if it is far enough from the others
end
l2 = setboundaries(l, [0, 5] => true, [5, 0] => true)  # Periodic boundary conditions
bonds = SiteDistance(<(1), l2)  # Connect sites that are closer than 1
plot(bonds)
plot!(l2)       # Unlike Bravais lattices, GenericLattice does not have default bonds
```

This tool is not very mature yet, but it can be useful for tasks not related to Bravais lattices. Note that boundary conditions or default nearest-neighbor hoppings can be set to the `GenericLattice` as well.