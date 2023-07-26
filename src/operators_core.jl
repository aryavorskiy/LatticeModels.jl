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
internal_basis(::LatticeBasis) = 1
internal_basis(b::CompositeLatticeBasis) = b.bases[1]
internal_basis(b::Basis) = throw(MethodError(internal_basis, (b,)))
internal_basis(any) = internal_basis(basis(any))

"""
    adjacency_matrix(op::Operator)

Generates an `AdjacencyMatrix` for the provided operator.
"""
function adjacency_matrix(op::LatticeOperator)
    matrix = Bool[!iszero(op.data[i, j])
                  for i in 1:length(lattice(op)), j in 1:length(lattice(op))]
    return AdjacencyMatrix(lattice(op), matrix)
end
function adjacency_matrix(op::CompositeLatticeOperator)
    n = length(internal_basis(op))
    ind(k) = (k - 1) * n + 1 : k * n
    matrix = Bool[!iszero(op.data[ind(i), ind(j)])
                  for i in 1:length(lattice(op)), j in 1:length(lattice(op))]
    return AdjacencyMatrix(lattice(op), matrix)
end

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

function tightbinding_hamiltonian(sample::Sample; tn=1, tnn=0, tnnn=0,
    field=NoField(), boundaries=BoundaryConditions())
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    internal_eye = one(internal).data
    for bond in default_bonds(l)
        add_hoppings!(builder, nothing, l, tn * internal_eye, bond, field, boundaries)
    end
    if tnn != 0
        for bond in default_nnbonds(l)
            add_hoppings!(builder, nothing, l, tnn * internal_eye, bond, field, boundaries)
        end
    end
    if tnnn != 0
        for bond in default_nnnbonds(l)
            add_hoppings!(builder, nothing, l, tnnn * internal_eye, bond, field, boundaries)
        end
    end
    return manybodyoperator(sample, to_matrix(builder))
end
tightbinding_hamiltonian(l::Lattice, b::Basis=GenericBasis(1); kw...) =
    tightbinding_hamiltonian(l ⊗ b; kw...)
