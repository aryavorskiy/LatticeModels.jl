We will begin the tutorial by defining the quantum system we are working with. This package introduces several abstractions to describe various quantum lattices. Let us observe them.

## Creating a lattice

The `Lattice` type is the base type for defining Bravais lattices. A Bravais lattice is by definition one or more sites forming a *unit cell* translated periodically by a set of vectors. The default lattice geometry in this package is a Bravais lattice with finite number of translations. This type of lattices will be further referred to as a *macrocell*.

Let's take a look at an example:

```@example env
using LatticeModels, Plots
l1 = SquareLattice(10, 5)                           # A 10x5 square lattice
l2 = HoneycombLattice(10, 5)                        # A 10x5 honeycomb lattice macrocell

p = plot(layout=2, size=(800, 350))                 # Create a plot with 2 subplots
plot!(p[1], l1, :pretty, title = "Square lattice")  # 'Pretty' lattice plot
plot!(p[2], l2, title = "Honeycomb lattice")        # Simple lattice plot
```

The positional arguments in the lattice constructors stand for number of translations along every lattice axis.
Note that `SquareLattice(5, 5, 5)` will generate a `5x5x5` 3D lattice, but `HoneycombLattice(5, 5, 5)` will throw an error, because a honeycomb lattice is 2D by its definition.

To achieve less trivial lattice geometry, you can exclude some sites from the lattice:
```@example env
macrocell = SquareLattice(10, 10)
center = macrocell[x = 5..6, y = 5..6]  # The four sites in the center

# All basic set operations work with `Lattice`s (see also `union`, `intersect`)
l = setdiff(macrocell, center) 
plot(l, :pretty)
```

Let's take a closer look at how we defined the center of the `macrocell`. `x = 5..6, y = 5..6` means that all 
sites of the `macrocell` with x-coordinate between 5 and 6 and y-coordinate also between 5 and 6 will get to the
`center`. This syntax allows selecting sites based on a set of parameters using these kwargs-style accessors:

- `x1`, `x2`, `x3` etc. stand for spatial coordinates. `x`, `y`, `z` are aliases for `x1`, `x2`, `x3`.
- `j1`, `j2`, `j3` etc. stand for *lattice axes*, in other words - the unit cell indices.
- `index` stands for number of the site in the unit cell.

Also note that instead of the `5..6` interval you can use any collection. Let's use all these techniques in one example:

```@example env
macrocell = HoneycombLattice(10, 10)
l1 = macrocell[j1 = 1:2:9, j2 = 2:2:10, index = 1]
l2 = macrocell[j1 = 2:2:10, j2 = 1:2:9, index = 2]

# Firstly, let's plot the macrocell with translucent markers
plot(macrocell, α = 0.3, lab = "macrocell")
# Then plot our resulting lattice
plot!(l1 ∪ l2, lab = "our lattice")
```

In some cases more complex lattice geometry is required. In such case we will need a different approach. But before that, let's get to know how lattice sites are treated within this package.

### Lattice sites

A [`LatticeSite`](@ref) is a struct describing where a site of some `Lattice` is located.
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
plot(macrocell, α = 0.3, lab = "macrocell")
plot!(l12, lab = "our lattice")
```

The function written using the do-syntax must return `true` if the site must be contained in the resulting lattice.
This syntax is more powerful, because it allows writing arbitrary code in the function body, which might be crucial in non-trivial cases.

## The `Sample`

The `Lattice` does not contain full information about the physical system - it is just a set of sites. 
Sometimes the system has additional degrees of freedom (for example, the spin of the particle).
Sometimes it has non-trivial boundary conditions. We need a separate type to store this information.

That is where the `Sample` comes in. If the boundaries are open, the `Sample` can be constructed by simple tensor product:

```julia
l = SquareLattice(10, 10)
spin = SpinBasis(1//2)
s = l ⊗ spin
```

Here `spin` is a `QuantumOptics.Basis` which describes the on-site phase state. It can be any other `QuantumOptics.Basis` - 
this package is built on top of `QuantumOpticsBase` and is fully compatible with its type ecosystem.

### Boundary conditions

All boundary conditions in this package are defined in terms of the macrocell: this means that if $R_i$ are the vectors
forming the supercell, boundary conditions can be defined in terms of $\psi(r), \psi(r + R_i)$. 
Three types of boundary conditions are supported:

- `PeriodicBoundary(axis)`: simple condition $\psi(r + R_{axis}) = \psi(r)$
- `TwistedBoundary(axis, Θ)`: condition $\psi(r + R_{axis}) = e^{i \theta} \psi(r)$
- `FunctionBoundary(f, axis)`: condition $\psi(r + R_{axis}) = f(r) \cdot \psi(r)$

Each condition affects only one axis. To use several conditions at once, use `BoundaryConditions`:

```julia
xbound = FunctionBoundary(1) do site
    # This construct may be useful in presence of magnetic field
    return exp(im * site.x) # Return the factor
end
ybound = TwistedBoundary(2, -pi / 2)
zbound = PeriodicBoundary(3)
bc = BoundaryConditions(xbound, ybound, zbound)
```

You can use a shorthand notation for boundary conditions like this:

```julia
f(site) = exp(im * site.x)
bc2 = BoundaryConditions(1 => f, 2 => -pi / 2, 3 => true)
```

After all this, you can add the boundary conditions to the sample:

```julia
s = Sample(l, spin, boundaries=bc)
```

### The `System`

The `Sample` contains no information about the particles - the quantity, the statistics etc. To do this, you need 
to create a `System`. It is simple:

```julia
sys1 = System(s, N = 3, statistics = BoseEinstein)  # Three bosons on the sample `s`
sys2 = System(s, statistics = FermiDirac, μ = 0)    # Non-interacting fermions with zero chemical potential
```

!!! tip
    Most functions in this package that create various operators accept a `Sample` or `System` as the first argument.
    You can avoid constructing it explicitly, because all these functions accept the lattice, internal basis and boundary conditions as separate arguments.

    The two function calls in this example will do the same thing:
    
    ```julia
    l = SquareLattice(10, 10)
    spin = SpinBasis(1//2)
    bs = BoundaryConditions(1 => true, 2 => true)
    sys = System(l, spin, boundaries = bs, μ = 0, statistics = FermiDirac)

    H1 = tightbinding_hamiltonian(sys)
    H2 = tightbinding_hamiltonian(l, spin, boundaries = bs)
    ```
