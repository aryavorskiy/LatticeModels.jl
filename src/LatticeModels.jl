module LatticeModels

using Reexport
@reexport using QuantumOpticsBase
@reexport using IntervalSets

@static if VERSION < v"1.8"
    allequal(seq) = all(==(first(seq)), seq)
end

const Nullable{T} = Union{Nothing,T}

include("core/sarrayutils.jl")
include("core/lattice.jl")
export lattice, sublattice, site_index, Coord
include("core/bonds.jl")
export adjacency_matrix, site_distance, SiteDistance, Translation
include("core/sites.jl")
include("core/latticevalue.jl")
export LatticeValue, siteproperty_value, coord_value, coord_values, project
include("core/recipes.jl")

include("lattices/bravais/unitcell.jl")
include("lattices/bravais/boundaries.jl")
export PeriodicBoundary, TwistedBoundary, FunctionBoundary, BoundaryConditions
include("lattices/bravais/lattice.jl")
export add_boundaries, UnitcellAxis, UnitcellIndex
include("lattices/bravais/bonds.jl")
export BravaisTranslation, NearestNeighbor
include("lattices/bravais/recipes.jl")

include("operators/bases.jl")
export LatticeBasis, ketstate, brastate
include("operators/system.jl")
export System, NParticles, FermiDirac, BoseEinstein, Hamiltonian
include("operators/magneticfield.jl")
export GaugeField, LineIntegralGaugeField
include("operators/builder.jl")
export OperatorBuilder
include("operators/constructoperator.jl")
export tightbinding_hamiltonian, construct_operator, construct_hamiltonian
include("operators/manybody.jl")
export interaction
include("operators/miscoperators.jl")
export siteproperty_operator, coord_operator, coord_operators
include("operators/latticeutils.jl")
export lattice_density, diag_reduce, apply_field!

include("spectrum.jl")
export Eigensystem, diagonalize, projector, densitymatrix, dos, ldos
include("currents.jl")
export Currents, currents_from, currents_from_to, mapgroup_currents
include("timesequence.jl")
export integrate, integrate!, differentiate, differentiate!, timestamps,
    TimeSequence
include("evolution.jl")
export @evolution

include("zoo/lattices.jl")
export SquareLattice, TriangularLattice, HoneycombLattice, KagomeLattice
include("zoo/magneticfields.jl")
export LandauGauge, SymmetricGauge, PointFlux
include("zoo/models.jl")
export hubbard, bosehubbard, fermihubbard, qwz, haldane, kanemele
include("zoo/currents.jl")
export DensityCurrents, OperatorCurrents

using Logging
try
    include("precompile.jl")
    _precompile_()
catch error
    @warn "Failed to precompile package due to unhandled exception:" error
end

end # module LatticeModels
