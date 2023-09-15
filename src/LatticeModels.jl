module LatticeModels

using Reexport
@reexport using QuantumOpticsBase
@reexport using IntervalSets

@static if VERSION < v"1.8"
    allequal(seq) = all(==(first(seq)), seq)
end

include("core/utils.jl")
include("core/abstract_lattice.jl")
include("core/adjacency.jl")
include("core/lattice_value.jl")
export LatticeValue, param_value, coord_values, project
include("core/plot_recipes.jl")

include("core/lattice_site.jl")
include("core/boundaries.jl")
export PeriodicBoundary, TwistedBoundary, FunctionBoundary, BoundaryConditions,
    PeriodicBoundaryConditions
include("core/lattice.jl")
export lattice, sublattice, site_index, site_distance
include("core/bonds.jl")
export SiteOffset, Bonds
include("core/bravais_plot_recipes.jl")

include("core/bases.jl")
export LatticeBasis, @increment, ketstate, brastate

include("core/field.jl")
export @field_def, NoField, MagneticField

include("core/sample.jl")
export  Sample, System, FermiDirac, BoseEinstein

include("operator_builder.jl")
export OperatorBuilder, Hamiltonian
include("operators_build.jl")
export tightbinding_hamiltonian, build_operator, build_hamiltonian
include("operators_manybody.jl")
export interaction
include("operators_utils.jl")
export coord_operator, coord_operators, lattice_density, diag_reduce, adjacency_matrix, apply_field!

include("spectrum.jl")
export Eigensystem, diagonalize, projector, densitymatrix, dos, ldos

include("currents.jl")
export materialize, currents_from, currents_from_to, map_currents

include("time_sequence.jl")
export init_record, integrate, integrate!, differentiate, differentiate!, timestamps,
    TimeSequence

include("evolution.jl")
    export @evolution

include("zoo.jl")
export SquareLattice, HoneycombLattice,
    LandauField, SymmetricField, FluxField,
    qwz, haldane, kanemele,
    DensityCurrents

using Logging
try
    include("precompile.jl")
    _precompile_()
catch error
    @warn "Failed to precompile package due to unhandled exception:" error
end

end # module LatticeModels
