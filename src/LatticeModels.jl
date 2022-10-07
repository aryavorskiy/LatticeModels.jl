module LatticeModels

include("lattice.jl")
export @lattice_def, LatticeIndex, Lattice, sublattice
export SquareLattice, HoneycombLattice

include("lattice_value.jl")
include("lattice_operator.jl")
export LatticeValue, LatticeVecOrMat, Basis,
    convert_inner_type, diag_operator, coord_operators, diag_aggregate, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!

include("hoppings.jl")
export hopping, hopping_operator, bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, projector, filled_projector

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, ChargeCurrents,
    materialize, current_lambda, lattice

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
