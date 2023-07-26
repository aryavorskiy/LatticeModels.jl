module LatticeModels

using Reexport
@reexport using QuantumOpticsBase
@reexport using IntervalSets

include("core/utils.jl")
include("core/lattice_site.jl")
include("core/lattice.jl")
export Lattice, sublattice, site_index, site_distance

include("core/lattice_value.jl")
export LatticeValue, coord_values, project

include("core/field.jl")
export @field_def, AbstractField, NoField

include("core/sample.jl")
export PeriodicBoundary, TwistedBoundary, FunctionBoundary, BoundaryConditions, Sample,
    FermiDirac, BoseEinstein

include("core/bonds.jl")
export Bonds

include("core/adjacency.jl")
export Domains, PairLhsGraph, PairRhsGraph

include("operators_core.jl")
export LatticeBasis, adjacency_matrix, tightbinding_hamiltonian
include("operators_build.jl")
export hoppings, build_hamiltonian
include("operators_manybody.jl")
export interaction
include("operators_utils.jl")
export coord_operators, coord, site_density

include("spectrum.jl")
export Eigensystem, diagonalize, projector, apply_to_eigenvalues, densitymatrix, dos, ldos

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, materialize, currents_from, currents_from_to, pairs_by_distance, map_currents

include("time_sequence.jl")
export init_record, integrate, integrate!, differentiate, differentiate!, time_domain,
    TimeSequence, LatticeValueSequence, LatticeArraySequence, CurrentsSequence

include("zoo.jl")
export SquareLattice, HoneycombLattice,
    LandauField, SymmetricField, FluxField,
    qwz, haldane,
    DensityCurrents

include("plot_recipes.jl")

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
