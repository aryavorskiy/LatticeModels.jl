module LatticeModels

using Reexport
@reexport using QuantumOpticsBase
@reexport using IntervalSets

@static if VERSION < v"1.8"
    allequal(seq) = all(==(first(seq)), seq)
end

include("core/sarrayutils.jl")
include("core/lattice.jl")
include("core/latticevalue.jl")
export LatticeValue, param_value, coord_values, project, @p_str
include("core/adjacency.jl")
include("core/recipes.jl")

include("lattices/bravais/site.jl")
include("lattices/bravais/boundaries.jl")
export PeriodicBoundary, TwistedBoundary, FunctionBoundary, BoundaryConditions,
    PeriodicBoundaryConditions
include("lattices/bravais/lattice.jl")
export lattice, sublattice, site_index, site_distance
include("lattices/bravais/bonds.jl")
export SiteOffset, Bonds
include("lattices/bravais/recipes.jl")

include("operators/bases.jl")
export LatticeBasis, @increment, ketstate, brastate
include("operators/system.jl")
export  Sample, System, FermiDirac, BoseEinstein
include("operators/magneticfield.jl")
export NoField, MagneticField
include("operators/constructor.jl")
export OperatorBuilder, Hamiltonian
include("operators/buildoperator.jl")
export tightbinding_hamiltonian, build_operator, build_hamiltonian
include("operators/manybody.jl")
export interaction
include("operators/latticedensity.jl")
export coord_operator, coord_operators, lattice_density, diag_reduce, adjacency_matrix

include("spectrum.jl")
export Eigensystem, diagonalize, projector, densitymatrix, dos, ldos
include("currents.jl")
export materialize, currents_from, currents_from_to, map_currents
include("timesequence.jl")
export init_record, integrate, integrate!, differentiate, differentiate!, timestamps,
    TimeSequence
include("evolution.jl")
export @evolution

include("zoo/lattices.jl")
export SquareLattice, HoneycombLattice
include("zoo/magneticfields.jl")
export LandauField, SymmetricField, FluxField
include("zoo/models.jl")
export qwz, haldane, kanemele
include("zoo/currents.jl")
export DensityCurrents

using Logging
try
    include("precompile.jl")
    _precompile_()
catch error
    @warn "Failed to precompile package due to unhandled exception:" error
end

end # module LatticeModels
