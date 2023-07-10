module LatticeModels

using Reexport
@reexport using QuantumOpticsBase

include("core/utils.jl")
include("core/lattice.jl")
export Bravais, LatticeSite, Lattice, sublattice, site_index, site_distance

include("core/lattice_value.jl")
export LatticeValue, coord_values, project

include("core/field.jl")
export @field_def, AbstractField, NoField

include("core/sample.jl")
export PeriodicBoundary, TwistedBoundary, FunctionBoundary, BoundaryConditions

include("core/bonds.jl")
export Bonds, Domains, PairLhsGraph, PairRhsGraph

include("operators.jl")
export LatticeBasis, diag_operator, coord_operators, apply_field!, coord, site_density, hoppings

include("spectrum.jl")
export Spectrum, spectrum, eigvals, eigvecs, projector, filled_projector,
    fermi_dirac, bose_einstein, dos, ldos

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
    TightBinding, SpinTightBinding, Haldane,
    DensityCurrents

include("plot_recipes.jl")

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
