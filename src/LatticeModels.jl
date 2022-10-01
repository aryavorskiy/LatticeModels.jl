module LatticeModels

include("lattice.jl")
export @lattice_def, LatticeIndex, AbstractLattice, SubLattice

include("lattice_value.jl")
include("lattice_operator.jl")
export LatticeValue, LatticeVecOrMat, Basis,
    convert_inner_type, diag_operator, coord_operators, diag_aggregate, @on_lattice

include("field.jl")
export @field_def, AbstractField

include("hoppings.jl")
export Hopping, BondSet, hopping_operator, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, projector, filled_projector

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, ChargeCurrents, MaterializedCurrents,
    materialize, current_lambda, lattice

end # module LatticeModels
