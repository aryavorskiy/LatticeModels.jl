module LatticeModels

include("lattice.jl")
export Bravais, LatticeSite, Lattice, dims, sublattice, site_coords, site_index,
    SquareLattice, HoneycombLattice

include("lattice_value.jl")
export LatticeValue, coord_values, project
include("lattice_operator.jl")
export LatticeOperator, Basis, ⊗, lattice, basis, dims_internal,
    diag_operator, coord_operators, diag_aggregate, ptrace, site_density, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!
export LandauField, SymmetricField, FluxField, NoField

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
export AbstractCurrents, DensityCurrents,
    materialize, current_lambda, lattice, pairs_by_distance, map_currents

include("precompile.jl")
_precompile_()

end # module LatticeModels
