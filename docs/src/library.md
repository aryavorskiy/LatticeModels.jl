# Internals

```@index
Modules = [LatticeModels]
```

## Lattice basics

```@autodocs
Modules = [LatticeModels]
Pages = ["core/lattice.jl", "lattices/genericlattice.jl", "lattices/bravais/unitcell.jl", 
    "lattices/bravais/lattice.jl"]
Filter = t -> t ∉ (LatticeModels.AbstractLattice, LatticeModels.AbstractSite)
```

## Lattice constructors

```@autodocs
Modules = [LatticeModels]
Pages = ["zoo/lattices.jl", "zoo/shapes.jl"]
```

## Bonds

```@autodocs
Modules = [LatticeModels]
Pages = ["core/bonds.jl", "lattices/bravais/bonds.jl", "lattices/bravais/nearestneighbor.jl"]
Filter = t -> t ∉ (LatticeModels.AbstractBonds, LatticeModels.AbstractTranslation, 
    LatticeModels.DirectedBonds)
```

## Boundary conditions

```@autodocs
Modules = [LatticeModels]
Pages = ["core/boundaries.jl"]
```

## LatticeValue

```@autodocs
Modules = [LatticeModels]
Pages = ["core/latticevalue.jl"]
```
 
## Operators and observables

```@autodocs
Modules = [LatticeModels]
Pages = ["operators/bases.jl", "operators/miscoperators.jl", "operators/manybody.jl", 
    "operators/latticeutils.jl"]
```

## Hamiltonians

```@autodocs
Modules = [LatticeModels]
Pages = ["operators/system.jl", "operators/builder.jl", "operators/constructoperator.jl",
    "operators/magneticfield.jl", "zoo/magneticfields.jl"]
```

## Built-in models
```@autodocs
Modules = [LatticeModels]
Pages = ["zoo/models.jl"]
```

## Diagonalization

```@autodocs
Modules = [LatticeModels]
Pages = ["spectrum.jl"]
```

## Green's function

```@autodocs
Modules = [LatticeModels]
Pages = ["greenfunction.jl"]
```

## Currents

```@autodocs
Modules = [LatticeModels]
Pages = ["currents.jl", "zoo/currents.jl"]
Filter = t -> t ∉ (LatticeModels.AbstractCurrents, LatticeModels.AbstractTranslation)
```

## Evolution

```@autodocs
Modules = [LatticeModels]
Pages = ["evolution.jl", "timesequence.jl"]
Filter = t -> t ∉ (LatticeModels.SchroedingerSolver,)
```

## Internals

```@docs
LatticeModels.AbstractLattice
LatticeModels.AbstractSite
LatticeModels.AbstractBonds
LatticeModels.DirectedBonds
LatticeModels.AbstractTranslation
LatticeModels.AbstractCurrents
LatticeModels.LookupTable
LatticeModels.SchroedingerSolver
LatticeModels.addlookuptable
```