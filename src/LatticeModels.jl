module LatticeModels

include("lattice.jl")
export Bravais, LatticeSite, Lattice, sublattice
export SquareLattice, HoneycombLattice

include("lattice_value.jl")
export LatticeValue, coord_values
include("lattice_operator.jl")
export LatticeOperator, Basis, ⊗, lattice, basis,
    diag_operator, coord_operators, diag_aggregate, ptrace, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!
export LandauField, SymmetricField, FluxField, NoField

include("hoppings.jl")
export hopping, hopping_operator, bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, eigvals, eigvecs, projector, filled_projector

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, DensityCurrents,
    materialize, current_lambda, lattice, is_near

# include("precompile.jl")
# _precompile_()

end # module LatticeModels
