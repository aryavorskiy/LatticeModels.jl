# Working with `LatticeValue`s

This chapter introduces the `LatticeValue` type, which describes a value defined on the sites of a lattice. `LatticeValue` implements the `AbstractArray` interface and is used to represent quantities such as the magnetization or local density of a state.

## Basics

The `LatticeValue` type is actually a wrapper around a `Vector` of values, where each value is associated with a site of the lattice. Working with `LatticeValue` is very similar to working with a `Vector`, but with some additional functionality.

```@example 1
using LatticeModels, Plots

l = SquareLattice(-2:2, -2:2)
x = coordvalue(l, :x)       # get the x coordinate of the lattice sites
v = zeros(Int, l)           # create a lattice value with zeros
v[x = 1 .. 2, y = 1 .. 2] = x   # set the value to the x coordinate on the top right
v[x = -2 .. 0] .= 3           # set the value to 3 on the left half of the lattice
r = randn(l)                # create a random lattice value (normal distribution)
v2 = v .^ 2 .+ r            # broadcast operations work as expected
heatmap(v2)
```

The `LatticeValue` type implements basic vector creation operations: `zero`, `zeros`, `one`, `ones`, `rand`, `randn`, `fill`, `copy`. Slices and views are also supported.

There are also other `AbstractVector` methods you can use:

```@repl 1
sum(abs2, v2)
extrema(v2)
argmax(v2)
ms = findall(x -> x < 3, v2)
v2[first(ms)]
```

### Indexing and slicing

Let's talk a bit more about the indexing and slicing of `LatticeValue`. There are several ways to index a `LatticeValue`. The return value in this case is either a scalar or another `LatticeValue` with a narrowed domain. The following indexing methods are supported:

- `v[site]` returns the value at site `site`.
- `v[[site1, site2, ...]]` returns a `LatticeValue` with the values at the specified sites. `site1`, `site2`, etc. are single sites grouped into an abstract array. The return value is a `LatticeValue` with the same values but narrowed to the specified sites.
- `v[lat]` returns a `LatticeValue` with the same values but narrowed to the sites of lattice `lat`. `lat` here must be a lattice, which is a subset of the lattice of `v`.
- `v[mask]` returns a `LatticeValue` with the same values but narrowed to the sites where `mask` is `true`. `mask` here must be a `LatticeValue` of `Bool` type, defined on the same (or a superset of the) lattice.
- `v[x = 1 .. 2, y = 1 .. 2]` returns a `LatticeValue` with the same values but narrowed to the sites where the x coordinate is in the range `1 .. 2` and the y coordinate is in the range `1 .. 2`. The keyword arguments here must be the names of the site parameters (see [Sites](@ref)), and the values can be any containers or single values. Pair notation is also supported: `v[Coord(1) => 1 .. 2, Coord(2) => 1 .. 2]`.

This indexing is valid for both reading and writing. Remember, however, that the right-hand side of the assignment must be a `LatticeValue` or a scalar value (in which case the destination site must be a single site, otherwise an error will be thrown). 

Therefore, `l[x=1, y=1] = 1` is not a valid assignment, because the left-hand side can contain multiple sites[^1] (for example, if it is a 3D lattice). But `l[!, x=1, y=1] = 1` is a valid assignment — adding `!` to the index means that the left-hand side is a single site.

[^1]: This is somewhat similar to the behavior of arrays in Julia: `v[1:1] = 1` will throw an error, even though the left-hand side is a single-element array.

!!! note
    The same indexing methods can be used to slice a lattice, a [`GreenFunction`](@ref), a [`Currents`](@ref) object or a [`TimeSequence`](@ref).

### Iteration and broadcasting

You can consider a `LatticeValue` as a vector of values, with its indices being the sites of the lattice. Therefore, 
iterating over `v` will yield the values of the `LaticeValue`, and `eachindex(v)` will return the lattice it is defined on.

Broadcasting operations work as expected. For example, `v .+ 1` will add 1 to each value of `v`, and `v .+ r` will add the corresponding values of `v` and `r`. However, there are some limitations: you cannot broadcast a `LatticeValue` with anything other than a scalar or another `LatticeValue`. Also the lattices must be the same, otherwise an error will be thrown.

```@repl 2
using LatticeModels
l = SquareLattice(4, 4);
x, y = coordvalues(l)
v = zeros(l)
v[x .< y] = x               # This will work
v[y = 1] .= 1               # This will work
v[x .> y] .= y              # This will not work - RHS on a different lattice
v[x = 1] .= [1, 2, 3, 4]    # This will not work - RHS is a vector
```

## Common operations

There are several common use-cases for the `LatticeValue` type. We will discuss some of them here.

### External parameter of a system

In many cases, you need to define a parameter that depends on the site that is used in the Hamiltonian. For example, the on-site potential in the tight-binding model. You can use `LatticeValue` for this purpose.

```julia
using LatticeModels
l = SquareLattice(4, 4)
v = zeros(l)
v[y = 0 .. 2] .= 1.0                          # add a potential barrier
H = tightbinding_hamiltonian(l, v, t1=-1.0) # create a tight-binding Hamiltonian
```

Here $v$ is a `LatticeValue` that represents the on-site potential. In other models you can use it to represent the magnetic field, for example: see [`qwz`](@ref).

### Wavefunctions

In some cases you need a custom-defined wavefunction. You can use `LatticeValue` for this purpose.

```julia
using LatticeModels
l = SquareLattice(10, 10)
x, y = coordvalues(l)
spin = SpinBasis(1//2)                  # create a spin basis
gauss = @. exp(-0.05 * ((x - 5.5) ^ 2 + (y - 5.5) ^ 2))
wave = @. exp(im * (x + y))             # create a plane wave
ψ = basisstate(spin, 1) ⊗ (@. gauss .* wave) + 
    basisstate(spin, 2) ⊗ (@. gauss * conj(wave))
```

Here `ψ` is a `QuantumOptics.Ket` wavefunction. In this example it is a superposition of two states with opposite spins and different momenta.

We will discuss this theme in more detail in the [States and Operators] section.

### Processing data

Many observables like local density are returned as a `LatticeValue`. You can process it quite easily.

```@example 2
using LatticeModels, Statistics
l = HoneycombLattice(Hexagon(), sites=120)
H = tightbinding_hamiltonian(l, t1=-1.0)    # create a tight-binding Hamiltonian
dens = localdensity(groundstate(H))         # calculate the local density of the ground state
r = shaperadius(l, Hexagon())               # get the radius of the lattice
bulk = HoneycombLattice(Hexagon(r * 0.8))
edge = setdiff(l, bulk)
println("Average bulk density: ", round(mean(dens[bulk]), digits=6))
println("Average edge density: ", round(mean(dens[edge]), digits=6))
```

The average local density in the bulk is much higher than on the edge, as expected.

We will discuss this theme in more detail in the [Measurements] section.

## Visualization

One key feature of `LatticeValue` is that it can be visualized using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package. There are several ways to do it.

The "classical" way is just to use the `plot` function. The result will be a scatter plot of the lattice sites with the value of the `LatticeValue` as the color and size of the markers.

Let us continue with the previous example:

```@example 2
using Plots
plot(dens)
```

This is the default behavior. You can customize the plot as usual with Plots.jl — for example, you can change the colormap, the marker size, etc. You can also set `markerscale=false` to disable the scaling of the marker size by the value of the `LatticeValue`:

```@example 2
# Dark theme makes everything look cooler
plot(dens, markerscale=false, markersize=12, c=:inferno, 
    title="Local density of the ground state", bg=:black)
```

Another way to visualize a `LatticeValue` is to use the shape plot `seriestype=:shape`, or `heatmap` function. This will create a tile plot of the lattice sites with the value of the `LatticeValue` as the color of the tiles.

```@example 2
heatmap(dens, title="Local density of the ground state")
```

You can also use pass `shape=:circle` to create a scatter plot with large circles instead of the default markers. The difference here is that the size of the circles will scale with the plot, unlike the markers in the `scatter` plot. `markerscale` is also supported here, but by default it is set to `false`.

Let's showcase all of these options:

```@example 2
p = plot(size=(1000, 850), layout=(2, 2))
heatmap!(p[1], dens, title="Hexagons, no scale")
heatmap!(p[2], dens, markerscale=true, title="Hexagons, scale")
heatmap!(p[3], dens, shape=:circle, title="Circles, no scale")
heatmap!(p[4], dens, shape=:circle, markerscale=true, title="Circles, scale")
```

!!! tip
    The shape plot is slower than the scatter plot, because it creates a separate shape for each site. If you
    are creating an animation or a large plot, you may want to use the scatter plot with custom-shaped markers instead.

Another important use case is dimension reduction. [Before](@ref Multi-dimensional-lattices) we already discussed how to plot a 2D slice of a 3D lattice. Here is an example of plotting a 1D slice of a 2D `LatticeValue`:

```@example 2
p = plot(size=(1000, 500), layout=(1, 2))
plot!(p[1], dens, st=:shape)
plot!(p[1], lattice(dens[j2 = 0]), :high_contrast)
plot!(p[2], dens[j2 = 0], axes=:x)
```

By projecting `axes=:x` the selected values on `j2 = 0` (e.g. the horizontal line in the middle of the plot) are shown as a 1D plot. Also we have shown the exact line where we took the slice from by plotting the markers with the `:high_contrast` setting.