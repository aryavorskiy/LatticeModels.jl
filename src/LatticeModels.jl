module LatticeModels

using Reexport
@reexport using QuantumOpticsBase
@reexport using IntervalSets

@static if VERSION < v"1.8"
    allequal(seq) = all(==(first(seq)), seq)
end

const Nullable{T} = Union{Nothing,T}

include("core/utils.jl")
include("core/lattice.jl")
export lattice, site_index, Coord
include("core/bonds.jl")
export AdjacencyMatrix, adjacencymatrix, sitedistance, adjacentsites, SiteDistance, Translation
include("core/boundaries.jl")
export BoundaryConditions, PeriodicBoundary, TwistedBoundary, FunctionBoundary, setboundaries
include("core/lookuptable.jl")
export addlookuptable
include("core/latticevalue.jl")
export LatticeValue, coordvalue, coordvalues, project
include("core/recipes.jl")

include("lattices/genericlattice.jl")
export GenericLattice, GenericSite

include("lattices/bravais/unitcell.jl")
export UnitCell, LatticeCoord, BasisIndex
include("lattices/bravais/lattice.jl")
export span_unitcells
include("lattices/bravais/bonds.jl")
export BravaisTranslation, Bravais
include("lattices/bravais/nearestneighbor.jl")
export NearestNeighbor, getnnbonds

include("operators/bases.jl")
export LatticeBasis, ketstate, brastate
include("operators/system.jl")
export System, NParticles, FermiDirac, BoseEinstein, Hamiltonian
include("operators/magneticfield.jl")
export GaugeField, LineIntegralGaugeField
include("operators/builder.jl")
export OperatorBuilder, FastOperatorBuilder
include("operators/constructoperator.jl")
export tightbinding_hamiltonian, construct_operator, construct_hamiltonian
include("operators/manybody.jl")
export interaction
include("operators/miscoperators.jl")
export coordoperator, coordoperators
include("operators/latticeutils.jl")
export localdensity, diag_reduce, apply_field!

include("spectrum.jl")
export Eigensystem, diagonalize, projector, densitymatrix, groundstate, findgroundstate
include("greenfunction.jl")
export GreenFunction, greenfunction, diagonalelements, dos, ldos
include("currents.jl")
export Currents, currentsfrom, currentsfromto, mapgroup_currents
include("timesequence.jl")
export integrate, integrate!, differentiate, differentiate!, timestamps,
    TimeSequence
include("evolution.jl")
export Evolution, CachedExp, KrylovKitExp

include("zoo/lattices.jl")
export @bravaisdef
export Chain, SquareLattice, TriangularLattice, HoneycombLattice, KagomeLattice, GrapheneRibbon
include("zoo/shapes.jl")
export scalefactor, shaperadius, fillshapes, addshapes!, deleteshapes!, removedangling!,
    BallND, Circle, Ball, Polygon, Triangle, Square, Hexagon, SiteAt, Rectangle, Path
include("zoo/magneticfields.jl")
export LandauGauge, SymmetricGauge, PointFlux
include("zoo/models.jl")
export hubbard, bosehubbard, fermihubbard, qwz, haldane, kanemele
include("zoo/currents.jl")
export DensityCurrents, LocalOperatorCurrents

include("recipes.jl")

using Logging
try
    include("precompile.jl")
    _precompile_()
catch error
    @warn "Failed to precompile package due to unhandled exception:" error
end

end # module LatticeModels
