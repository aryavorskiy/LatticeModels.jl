This section contains the documentation for the internal structure of `LatticeModels.jl`. 

# Advanced features

These features can be useful in non-trivial cases, but are not necessary for basic usage.

## Diagonalizing Hamiltonians

## `SchroedingerSolver`s

## Currents

# Interfaces

The base of `LatticeModels.jl` is its interfaces, allowing to define lattices with arbitrary geometry, topology and boundary conditions. The interfaces are designed to be as flexible as possible, allowing to define new types of lattices and bonds.

## `AbstractLattice`

### Basic functions

### Mutable lattices

### Site parameters

### Lattice metadata

### Site lookup

### Shapes

## `AbstractBonds`

### Three main types of bonds

### Adapting bonds to the lattice

### Boundary conditions