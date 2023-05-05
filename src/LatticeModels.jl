module LatticeModels

include("utils.jl")
include("lattice.jl")
export Bravais, LatticeSite, Lattice, dims, sublattice, site_coords, site_index, site_distance

include("lattice_value.jl")
export LatticeValue, coord_values, project

include("lattice_array.jl")
include("lattice_basis.jl")
export LatticeArray, LatticeOperator, Basis, ⊗, lattice, basis, dims_internal,
    diag_operator, coord_operators, diag_reduce, ptrace, site_density, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!, NoField

include("hoppings.jl")
export hopping, hopping_operator,
    DomainsSelector, PairLhsSelector, PairRhsSelector,
    bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, eigvals, eigvecs, projector, filled_projector,
    fermi_dirac, bose_einstein, dos, ldos

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, materialize, current_lambda, lattice, pairs_by_distance, map_currents

include("record.jl")
export init_record, integrate, time_domain,
    LatticeRecord, LatticeValueRecord, LatticeArrayRecord, CurrentsRecord

include("zoo.jl")
export SquareLattice, HoneycombLattice,
    LandauField, SymmetricField, FluxField,
    TightBinding, SpinTightBinding, Haldane,
    DensityCurrents

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
