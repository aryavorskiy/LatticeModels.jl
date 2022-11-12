module LatticeModels

include("lattice.jl")
export Bravais, LatticeSite, Lattice, sublattice, site_coords, site_index,
    SquareLattice, HoneycombLattice

include("lattice_value.jl")
export LatticeValue, coord_values
include("lattice_operator.jl")
export LatticeOperator, Basis, ⊗, lattice, basis,
    diag_operator, coord_operators, diag_aggregate, ptrace, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!
export LandauField, SymmetricField, FluxField, NoField

include("hoppings.jl")
export hopping, hopping_operator, pairs_by_domains, pairs_by_lhs, pairs_by_rhs,
    bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, eigvals, eigvecs, projector, filled_projector

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, DensityCurrents,
    materialize, current_lambda, lattice, pairs_by_adjacent, pairs_by_distance

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
