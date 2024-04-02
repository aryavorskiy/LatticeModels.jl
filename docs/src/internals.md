This section contains the documentation for the internal structure of `LatticeModels.jl`. 

!!! warning
    This section of documentation is under construction. Some parts may be incomplete.

# Advanced features

These features can be useful in non-trivial cases, but are not necessary for basic usage.

## `AbstractLattice` interface

The base of `LatticeModels.jl` is its interfaces, allowing to define lattices with arbitrary geometry, topology and boundary conditions. The [`LatticeModels.AbstractLattice`](@ref) interface is the main interface for defining lattices.

Generally speaking, a lattice is a set of sites. Each site, in turn, has its spatial coordinates in its `coords` field and maybe some additional properties. It also must be a subtype of [`LatticeModels.AbstractSite`](@ref).

Note that the bonds between sites and the boundary conditions are initially not part of the lattice, but are added to its metadata later.

### Basic functions

### Site lookup

### Mutable lattices

### Site properties

### Lattice metadata

### Shapes

## `AbstractBonds` interface

The [`LatticeModels.AbstractBonds`](@ref) interface is used to define different types of bonds between sites. Most generally speaking,
such object is a mapping that decides if the sites are connected for each pair.

Note that there are three basic types of bonds in `LatticeModels.jl`:
- `LatticeModels.AbstractBonds`: a most general interface. Basically, it is just a mapping from site pairs to boolean values.
- [`LatticeModels.DirectedBonds`](@ref): this type of bonds defines a set of bonds that has a defined direction. The whole topology can be defined by the "destination" sites for each site. Since the bonds are usually sparse, the general performance of this type of bonds is much higher.
- [`LatticeModels.AbstractTranslation`](@ref): this is a subtype of `DirectedBonds`, where every site has one or zero "destination" sites. This allows to increase the performance even more, and also to transform the sites in a convenient manner:

```julia
site1 = lat[!, x = 1, y = 1]    # Get the site at [1, 1]
T = Translation(lat, [1, 0])    # Translate the site by [1, 0] vector
site2 = site1 + T               # `site2` is at [2, 1]
```

### Adapting bonds to the lattice

### Boundary conditions

## Diagonalizing the Hamiltonian

It is very easy to diagonalize a matrix in Julia. However, problems can arise when the matrix is of some custom type
(e. g. sparse or a GPU array). By default `LatticeModels.jl` makes use of `KrylovKit.jl` to solve the eigenproblem using the Lanczos algorithm for non-trivial matrix types. However, sometimes it is necessary to use a different algorithm. The `LatticeModels.diagonalize_routine` is a simple way to add a new algorithm to the default toolchain.

## `SchroedingerSolver`s

The [`LatticeModels.SchroedingerSolver`](@ref) interface is used to solve the time-dependent Schrödinger equation. It is used in the `Evolution` struct to perform unitary evolution. As with the diagonalization problem, one can add a new algorithm to the default toolchain by creating a new `SchroedingerSolver` type.

## Currents

The [`LatticeModels.AbstractCurrents`](@ref) interface allows to define different types of currents on the lattice. This allows it to be a lazy object, which computes the currents only when needed.

To implement basic currents semantics, you need to define the following methods:
- `LatticeModels.lattice(your_currents)`: returns the lattice, on which the currents are defined.
- `Base.getindex(your_currents, i::Int, j::Int)`: returns the current between sites with numbers `i` and `j`. This is done in such a manner, because you do not usually need the site properties to calculate the currents.