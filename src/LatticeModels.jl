module LatticeModels

include("lattice.jl")
export @lattice_def, LatticeIndex, AbstractLattice, SubLattice, sublattice

include("lattice_value.jl")
include("lattice_operator.jl")
export LatticeValue, LatticeVecOrMat, Basis,
    convert_inner_type, diag_operator, coord_operators, diag_aggregate, @on_lattice

include("field.jl")
export @field_def, AbstractField, apply_field!

include("hoppings.jl")
export Hopping, hopping_operator, BondSet, bonds, is_adjacent, @hopping_operator

include("hamiltonian.jl")
export @hamiltonian, Spectrum, spectrum, projector, filled_projector

include("evolution.jl")
export @evolution

include("currents.jl")
export AbstractCurrents, ChargeCurrents, MaterializedCurrents,
    materialize, current_lambda, lattice

macro precompile(expr)
    !Meta.isexpr(expr, :call) && error("function call expected")
    fn_name, fn_args... = expr.args
    call_str = "Failed to precompile " * string(expr)
    for i in eachindex(fn_args)
        !Meta.isexpr(fn_args[i], :(::), 1) && error("type definition (e. g. ::MyType) expected")
        fn_args[i] = esc(only(fn_args[i].args))
    end
    types_tpl = Expr(:tuple, fn_args...)
    :(!precompile($(esc(fn_name)), $types_tpl) && @warn($call_str))
    :()
end
POSSIBLE_LATTYPES = SquareLattice, SubLattice{SquareLattice}
# POSSIBLE_LATTYPES = HoneycombLattice, SubLattice{HoneycombLattice}
POSSIBLE_ELTYPES = Int, Float64, ComplexF64

# !!! Great precompile impact
_sf(site) = true
FT = typeof(_sf)
for LT in POSSIBLE_LATTYPES
    LVT = LatticeVecOrMat{LT, Matrix{ComplexF64}}
    for ET in POSSIBLE_ELTYPES
        HT = Hopping{Matrix{ET}}
        @precompile diag_operator(::LT, ::Matrix{ET})
        @precompile hopping_operator(::FT, ::LT, ::HT)
        @precompile hopping_operator(::LT, ::HT)
        @precompile projector(::Spectrum{LT, Matrix{ET}})
        for math_op in (*, /, ^)
            @precompile math_op(::LVT, ::ET)
        end
        @precompile bonds(::LT, ::HT)
        @precompile bonds(::LT, ::HT, ::Vararg{HT})
    end
    @precompile _propagate_lattice_args(::FT, ::LT)
    @precompile iterate(::LT)
    @precompile iterate(::LT, ::Tuple{LatticeIndex, Int, Int})
    @precompile sublattice(::FT, ::LT)
    @precompile coord_operators(::Basis{LT})
    @precompile diag_operator(::FT, ::LT)
    @precompile diag_aggregate(::FT, ::LVT)
    @precompile spectrum(::LVT)
    @precompile projector(::FT, ::LVT)
    @precompile bonds(::LVT)
    @precompile _unwrap(typeof(+), NTuple{2, LVT})
    @precompile _unwrap(typeof(*), Tuple{Float64, LVT})
    @precompile ^(::BondSet{LT}, ::Int)
    @precompile is_adjacent(::BondSet{LT})
    @precompile is_adjacent(::BondSet{LT}, ::LatticeIndex, ::LatticeIndex)
end

@precompile _hamiltonian_block(::Expr)
@precompile _evolution_block(::Expr, ::Expr)
@precompile _wrap_smart!(::Expr)

end # module LatticeModels
