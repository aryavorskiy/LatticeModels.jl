using SparseArrays
import QuantumOpticsBase
import QuantumOpticsBase: Basis, SparseOperator, SparseOpPureType, coefficient

struct LatticeBasis{LT<:Lattice} <: QuantumOpticsBase.Basis
    shape::Int
    latt::LT
    LatticeBasis(l::LT) where LT<:Lattice = new{LT}(length(l), l)
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt

onebodybasis(sample::Sample) = sample.internal ⊗ LatticeBasis(sample.latt)
onebodybasis(sample::LatticeSample) = LatticeBasis(sample.latt)

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::LatticeSite) =
    basisstate(T, b, site_index(b.latt, site))

const CompositeLatticeBasis{S, BT, LT} = CompositeBasis{S, Tuple{BT, LatticeBasis{LT}}}
const AbstractLatticeBasis = Union{LatticeBasis, CompositeLatticeBasis}

const LatticeOperator{MT} = Operator{BT, BT, MT} where BT<:LatticeBasis
const CompositeLatticeOperator{MT} = Operator{BT, BT, MT} where BT<:CompositeLatticeBasis
const AbstractLatticeOperator{MT} = Operator{BT, BT, MT} where BT<:AbstractLatticeBasis

lattice(lb::LatticeBasis) = lb.latt
lattice(b::CompositeLatticeBasis) = lattice(b.bases[2])
lattice(b::Basis) = throw(MethodError(lattice, (b,)))
lattice(any) = lattice(basis(any))
internal_basis(::LatticeBasis) = throw(ArgumentError("Lattice basis has no internal"))
internal_basis(b::CompositeLatticeBasis) = b.bases[1]
internal_basis(b::Basis) = throw(MethodError(internal_basis, (b,)))
internal_basis(any) = internal_basis(basis(any))
internal_length(any) = internal_length(basis(any))
internal_length(::LatticeBasis) = 1
internal_length(b::Basis) = length(internal_basis(b))

function add_diagonal!(builder, op, diag)
    for i in 1:length(diag)
        increment!(builder, op, i, i, diag[CartesianIndex(i)])
    end
end

@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
@inline _get_bool_value(g::AbstractGraph, ::Lattice, site1::LatticeSite, site2::LatticeSite) =
    match(g, site1, site2)

function add_hoppings!(builder, selector, l::Lattice, op, bond::Bonds,
        field::AbstractField, boundaries::BoundaryConditions)
    dims(bond) > dims(l) && error("Incompatible dims")
    trv = radius_vector(l, bond)
    for site1 in l
        lp = site1 + bond
        lp === nothing && continue
        add_hoppings!(builder, selector, l, op, site1 => LatticeSite(lp, site1.coords + trv), field, boundaries)
    end
end

function add_hoppings!(builder, selector, l::Lattice, op, bond::SingleBond,
        field::AbstractField, boundaries::BoundaryConditions)
    site1, site2 = bond
    p1 = site1.coords
    p2 = site2.coords
    factor, site2 = shift_site(boundaries, l, site2)
    i = @inline site_index(l, site1)
    j = @inline site_index(l, site2)
    i === nothing && return
    j === nothing && return
    !_get_bool_value(selector, l, site1, site2) && return
    total_factor = exp(-2π * im * line_integral(field, p1, p2)) * factor
    !isfinite(total_factor) && error("got NaN or Inf when finding the phase factor")
    increment!(builder, op, i, j, total_factor)
    increment!(builder, op', j, i, total_factor')
end
