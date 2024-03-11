We will begin the tutorial by defining the quantum system we are working with. This package introduces several abstractions to describe various quantum lattices. Let us observe them.

## Creating a lattice

The `BravaisLattice` type is the base type for defining Bravais lattices. A Bravais lattice is by definition one or more sites forming a *unit cell* translated periodically by a set of vectors. The default lattice geometry in this package is a Bravais lattice with finite number of translations.

Let's take a look at an example:

```@example env
using LatticeModels, Plots
l1 = SquareLattice(10, 5)                           # A 10x5 square lattice
l2 = HoneycombLattice(10, 5)                        # A 10x5 honeycomb lattice

p = plot(layout=2, size=(800, 350))                 # Create a plot with 2 subplots
plot!(p[1], l1, :pretty, title = "Square lattice")  # 'Pretty' lattice plot
plot!(p[2], l2, title = "Honeycomb lattice")        # Simple lattice plot
```

The positional arguments in the lattice constructors stand for number of translations along every lattice axis.
Note that `SquareLattice(5, 5, 5)` will generate a `5x5x5` 3D lattice, but `HoneycombLattice(5, 5, 5)` will throw an error, because a honeycomb lattice is 2D by its definition.

To achieve less trivial lattice geometry, you can exclude some sites from the lattice:
```@example env
latt = SquareLattice(10, 10)
center = latt[x = 5..6, y = 5..6]  # The four sites in the center

# All basic set operations work with `Lattice`s (see also `union`, `intersect`)
l = setdiff(latt, center) 
plot(l, :pretty)
```

Let's take a closer look at how we defined the center of the `latt`. `x = 5..6, y = 5..6` means that all 
sites of the `latt` with x-coordinate between 5 and 6 and y-coordinate also between 5 and 6 will get to the
`center`. This syntax allows selecting sites based on their coordinates. Note that instead of the `5..6` interval you can use any collection. 

You can also select sites by other parameters like unit cell indices using similar syntax:

```@example env
center2 = latt[p"x" => 5..6, p"y" => 5..6]  # The same as `center`

base_latt = HoneycombLattice(10, 10)
l1 = base_latt[p"j1" => 1:2:9, p"j2" => 2:2:10, p"index" => 1]
l2 = base_latt[p"j1" => 2:2:10, p"j2" => 1:2:9, p"index" => 2]

# Firstly, let's plot the `base_latt` with translucent markers
plot(base_latt, α = 0.3, lab = "base_latt")
# Then plot our resulting lattice
plot!(l1 ∪ l2, lab = "our lattice")
```

This syntax is called _site parameters_ and is generally more preferred, because it is infinitely flexible and can be redefined for any site type. A site parameter is a string literal preceded by `p` character:

- `p"x$i"` stands for `i`-th coordinate axis. `p"x"`, `p"y"`, `p"z"` are aliases for the first three coordinate axes.
- `p"j$i"` stands for `i`-th lattice axis (valid for Bravais sites only).
- `p"index"` stands for the site's index in the Bravais lattice basis.

!!! tip
    This syntax is used not only for constructing sublattices, but also:
    - Generating and accessing [`LatticeValue`s](@ref processing_results).
    - Creating and [diagonal operators](@ref diagonal_operators).
    - Building projections on different axes.
    
    Note that you can use `:x` instead of `p"x"`, `:y` instead of `p"y"` and `:z` instead of `p"z"` in any context.

In some cases more complex lattice geometry is required. In such case we will need a different approach. But before that, let's get to know how lattice sites are treated within this package.

### Lattice sites

A [`BravaisSite`](@ref) is a struct describing where a site of some `BravaisLattice` is located.
It has three fields:
- `unit_cell`: the location of the unit cell containing the site
- `index`: number of the site in the unit cell
- `coords`: the spatial coordinates

You can consider a `Lattice` an array of `LatticeSite`s: accessing a `Lattice` at an integer index or iterating over it 
yields a `LatticeSite`. You can also obtain a site using the same kwargs-style accessors.

```julia
l = SquareLattice(5, 5)
site1 = l[7]            # Accessing a site by its number
site2 = l[x = 2, y = 2] # The same site by coordinates
site1 == site2          # true
```

!!! tip
    Plot your lattice with the `:pretty` modifier to see the number of each site.

Now when you know about how `LatticeSite`s work, the last example can be rewritten:

```@example env
l12 = HoneycombLattice(10, 10) do site
    j1, j2 = site.unit_cell
    if j1 % 2 == 1 && j2 % 2 == 0
        return site.index == 1
    elseif j1 % 2 == 0 && j2 % 2 == 1
        return site.index == 2
    else return false
    end
end
plot(base_latt, α = 0.3, lab = "base_latt")
plot!(l12, lab = "our lattice")
```

The function that uses the `do`-syntax must return `true` if the site must be contained in the resulting lattice.
This syntax is more powerful, because it allows writing arbitrary code in the function body, which might be crucial in non-trivial cases.

### Boundary conditions

All boundary conditions in this package are defined in terms of unit cell vectors: this means that if $R_i$ are the vectors
forming the unit cell, boundary conditions can be defined in terms of $t_j$ vector that connects the wavefunction in points $\psi(r), \psi(r + \sum_j t_j R_j)$. 
Three types of boundary conditions are supported:

- `PeriodicBoundary(t)`: simple condition $\psi(r + \sum_j t_j R_j) = \psi(r)$
- `TwistedBoundary(t, Θ)`: condition $\psi(r + \sum_j t_j R_j) = e^{i \theta} \psi(r)$
- `FunctionBoundary(f, t)`: condition $\psi(r + \sum_j t_j R_j) = f(r) \cdot \psi(r)$

Each condition defines one translation axis. To define several at once, use `BoundaryConditions`:

```julia
xbound = FunctionBoundary([10, 0, 0]) do site
    # This construct may be useful in presence of magnetic field
    return exp(im * site.x) # Return the factor
end
ybound = TwistedBoundary([0, 10, 0], -pi / 2)
zbound = PeriodicBoundary([0, 0, 10])
bc = BoundaryConditions(xbound, ybound, zbound)
```

You can use a shorthand notation for boundary conditions like this:

```julia
f(site) = exp(im * site.x)
bc2 = BoundaryConditions([10, 0, 0] => f, [0, 10, 0] => -pi / 2, [0, 0, 10] => true)
```

After all this, you can add the boundary conditions to the lattice constructor:

```julia
l = SquareLattice(10, 10, 10, boundaries=bc2)
# Same shorthand applies
l2 = SquareLattice(10, 10, 10, boundaries=
    ([10, 0, 0] => f, [0, 10, 0] => -pi / 2, [0, 0, 10] => true))
```

### The `System`

The `Lattice` does not contain full information about the physical system - it is just a set of sites. 
Sometimes the system has additional degrees of freedom (for example, the spin of the particle).
Sometimes it has non-trivial boundary conditions. We need a separate type to store this information.

That is where the `System` comes in. If there is only one particle on the lattice, the `System` can be constructed by simple tensor product:

```julia
l = SquareLattice(10, 10)
spin = SpinBasis(1//2)
s = l ⊗ spin
```

Here `spin` is a `QuantumOptics.Basis` which describes the on-site phase state. It can be any other `QuantumOptics.Basis` - 
this package is built on top of `QuantumOpticsBase` and is fully compatible with its ecosystem.

You can build more complex `System`s over this one-particle abstraction:

```julia
sys1 = System(s, μ = 0, statistics = FermiDirac, T = 0.1)   # Non-interacting fermions with zero chemical potential
sys2 = System(s, N = 3, statistics = BoseEinstein, T = 0)   # Three non-interacting bosons on the sample `s`
sys3 = NParticles(s, 3, statistics = BoseEinstein, T = 1)   # Three interacting bosons on the sample `s`
```

The temperature can be set using the `T` keyword - it will be used to calculate the density matrix of the system.

!!! tip
    Note that you can use `mu` keyword instead of `μ` if your environment does not support unicode input.

!!! tip
    Most functions in this package that create various operators accept a `System` as the first argument.
    For non-interacting systems you can avoid constructing it explicitly, because all these functions accept the lattice, internal basis, boundary conditions etc. as separate arguments. 
    You can also pass boundary conditions to the `System` constructor to override existing boundary conditions of `l` lattice.

    The two `tightbinding_hamiltonian` function calls in this example will do the same thing:
    
    ```julia
    l = SquareLattice(10, 10)
    spin = SpinBasis(1//2)
    bs = BoundaryConditions([10, 0] => true, [0, 10] => true)
    sys = System(l, spin, boundaries = bs, μ = 0, statistics = FermiDirac)

    H1 = tightbinding_hamiltonian(sys)
    H2 = tightbinding_hamiltonian(l, spin, boundaries = bs, μ = 0, statistics = FermiDirac)
    ```