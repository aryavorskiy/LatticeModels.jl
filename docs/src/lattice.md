## Simple Bravais lattice

The simplest variant of a finite Bravais lattice is a *macro cell*, 
which is the lattice basis translated finite number of times along every translation vector.

```@setup env
using LatticeModels, Plots, LinearAlgebra
```

Constructing a macro cell is simple: the lattice constructor accepts positional arguments
which will be translation ranges along all dimensions.

Note that lattices of some types can be of any dimensionality, while others can not.

```jldoctest; setup=:(using LatticeModels)
julia> SquareLattice(10, 10)
100-site square lattice on 10×10 macro cell

julia> SquareLattice(3, 3, 3)
27-site square lattice on 3×3×3 macro cell

julia> HoneycombLattice(5, 5)
50-site honeycomb lattice on 5×5 macro cell (2-site basis)

julia> HoneycombLattice(3, 3, 2)
ERROR: MethodError: no method matching HoneycombLattice(::Int64, ::Int64, ::Int64)
[...]
```

A lattice can be scatter-plotted to see how its sites are located and which index is assigned to each site:

```@example env
p = plot(size=(800, 350), layout=2)
plot!(p[1], SquareLattice(10, 5))
plot!(p[2], HoneycombLattice(8, 4))
```

## [Lattice sites and axis descriptors](@id axis_descriptors)

A [`LatticeSite`](@ref) is a struct describing where a site of some Bravais lattice is located: 
it stores the location of the unit cell and the site's number in the lattice basis.
No information about its spatial coordinates is stored explicitly, but it can be obtained by passing the corresponding `Lattice` and the site to the [`site_coords`](@ref) function.

Iterating over any `Lattice` will yield `LatticeSite`s. 
You can also get them by indexing the `Lattice` object with integers or coordinate keywords:

```@example env
l = SquareLattice(5, 5)
site1 = l[7]
site2 = l[x = 2, y = 2]
site1 == site2
```

These keywords are called *axis descriptors*. There are two types of them:
- *Coordinate axes* are described by integers or symbols of form `:x$i` (where `$i` is the integer index of the axis). For indices 1-3 descriptors `:x`, `:y` and `:z` respectively are also allowed. Use these to select lattice sites by coordinate.
- *Lattice axes* parallel to lattice translation vectors are described by symbols of form `:j$i` (where `$i` is the integer index of the axis). Use these to select lattice sites by unit cell index.
- A special *basis index descriptor* `:index` can be used to select sites by basis index.
Axis descriptors can also be used for selecting sites (or sublattices) by lattice indices. Take a look at the example:

```@example env
l = HoneycombLattice(5, 5)
plot(l, leg=:bottomright)
plot!(l[j1=2], pretty=false, lab="Select a unitcell row")
plot!(l[index=1], pretty=false, ms=10, alpha=0.3, lab="Select a Bravais sublattice")
```

## Sublattices

Suppose we want to create a lattice with non-trivial geometry (for example, with holes). 
This can be done by deleting some of the sites from the macro cell. There are three ways ways to do this:

**The Convenient way**

```@example env
l = SquareLattice(10, 10)

x, y = coord_values(l)
l1 = l[@. abs(x) > 3 || abs(y) > 3]
```

Here we first create the macro cell, then find the coordinate values for its sites.
After that we use LatticeValue broadcasting, see [Lattice values](@ref) for more detail.

!!! tip
    This way to define sublattices is preferred, because code like this is the most readable.
    The `x, y` coordinate values will be also helpful to create other sublattices or slices.

**The Low-level way**

```@example env
l2 = sublattice(l) do site, (x, y); abs(x) > 3 || abs(y) > 3; end
```

The lambda must accept a LatticeSite and a coordinate vector as positional arguments 
and return whether the site should be included or not.

**The In-place way**

```@example env
l3 = SquareLattice(10, 10) do site, (x, y); abs(x) > 3 || abs(y) > 3; end
```

This notation is exactly the same as the low-level way, but done in one line.

!!! warning
    This notation should be used only if you need to access site indices or if the `x, y` coordinate values will not be needed further in the program.

The plot recipe for sublattices shows excluded sites with translucent markers by default and also prints out integer indices for all included sites. Pass keyword argument `pretty=false` to suppress this behavior.

```@example env
l = HoneycombLattice(6, 6)
x, y = coord_values(l)
plot(l[@. -5 < x * √3 + y < 5.5])
```

## Custom lattice types

It is quite likely that you will need more types of lattices than this package provides by default. In such cases you need to define a new type. Follow these steps:

**Create an exact alias**

Select a `Symbol` that will be the LatticeSym for this type and define an alias for the `Lattice{LatticeSym, N, NB}` type. The alias must not have any type parameters except for the dimension count if needed:

```julia
const HoneycombLattice = Lattice{:honeycomb, 2, 2}
const SquareLattice{N} = Lattice{:square, N, 1}
```

**Define the constructor**

The only positional arguments allowed are the macro cell size.[^1] The [`Bravais`](@ref) object must be generated in the constructor and passed to the default constructor `Lattice(sym, sz, bvs)`.

[^1]: This is done with purpose to achieve code consistency. Also in-place sublattice generation will almost certainly be broken. Use keyword arguments if you need additional parameters for some lattice type.

Let us define our own lattice type:
```@example env
const GrapheneLattice = Lattice{:graphene, 3, 2}
function GrapheneLattice(sz::Vararg{Int, 3})
    bvs = Bravais([1 1/2 0; 0 √3/2 0; 0 0 2], [0 1/2; 0 √3/6; 0 0])
    Lattice(:graphene, sz, bvs)
end
nothing # hide
```

A sublattice constructor will be generated by default:

```@example env
gl = GrapheneLattice(6, 6, 3) do site, (x, y, z)
    -5 < x * √3 + y < 5.5
end
plot(gl, pretty=false)
```

!!! warning
    Please note that if the type alias is dimension-parametric, you must define the constructor *for a concrete type*, otherwise you will almost definitely break the lattice constructor dispatch:
    ```julia
    SquareLattice(sz::Vararg{Int, N}) where N = ...     # Wrong!
    SquareLattice{N}(sz::Vararg{Int, N}) where N = ...  # Correct
    ```
    In the second example `SquareLattice(sz::Vararg{Int, N}) where N` constructor will be generated automatically.