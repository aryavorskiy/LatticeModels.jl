module LatticeModels

using Reexport
@reexport using QuantumOpticsBase

include("utils.jl")
include("lattice.jl")
export Bravais, LatticeSite, Lattice, dims, sublattice, site_index, site_distance

include("lattice_value.jl")
export LatticeValue, coord_values, project


include("lattice_basis.jl")
export LatticeArray, LatticeOperator, Basis, ⊗, lattice, basis, dims_internal,
    diag_operator, coords, diag_reduce, ptrace, site_density, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!, NoField

include("hoppings.jl")
export hopping, hoppings,
    DomainsSelector, PairLhsSelector, PairRhsSelector,
    bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, eigvals, eigvecs, projector, filled_projector,
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

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
