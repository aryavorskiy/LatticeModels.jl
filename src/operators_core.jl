using SparseArrays
import QuantumOpticsBase

struct LatticeBasis{LT<:Lattice} <: QuantumOpticsBase.Basis
    shape::SVector{1, Int}
    latt::LT
    LatticeBasis(l::LT) where LT<:Lattice = new{LT}(SA[length(l)], l)
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt
lattice(lb::LatticeBasis) = lb.latt

function one_particle_basis(sample::Sample)
    lb = LatticeBasis(sample.latt)
    length(sample.internal) == 1 ? lb : sample.internal ⊗ lb
end

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::LatticeSite) =
    basisstate(T, b, site_index(b.latt, site))

const LatticeOperator{BT, MT} = Operator{BT, BT, MT} where BT<:LatticeBasis

struct SparseMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    SparseMatrixBuilder{T}(sz) where T = new{T}(sz, [], [], [])
    SparseMatrixBuilder{T}(sz...) where T = SparseMatrixBuilder{T}(sz)
end
Base.size(smb::SparseMatrixBuilder) = smb.size
to_matrix(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int, factor=1)
    push!(builder.Is, i1)
    push!(builder.Js, i2)
    push!(builder.Vs, rhs * factor)
end

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::AbstractMatrix, i1::Int, i2::Int, factor=1)
    N = size(rhs)[1]
    for i in 1:N, j in 1:N
        v = rhs[i, j]
        iszero(v) && continue
        increment!(builder, v, i + N * (i1 - 1), j + N * (i2 - 1), factor)
    end
end

function add_diagonal!(builder, op, diag)
    for i in 1:length(diag)
        increment!(builder, op, i, i, diag[CartesianIndex(i)])
    end
end

"""
    radius_vector(l::Lattice, hop::Hopping)
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macrocell to minimize the distance between them).
"""
function radius_vector(l::Lattice, hop::Bonds)
    i, j = hop.site_indices
    bravais(l).basis[:, j] - bravais(l).basis[:, i] +
     mm_assuming_zeros(bravais(l).translation_vectors, hop.translate_uc)
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
        site1.basis_index != bond.site_indices[1] && continue
        site2 = LatticeSite(add_assuming_zeros(site1.unit_cell, bond.translate_uc),
            bond.site_indices[2], site1.coords + trv)
        add_hoppings!(builder, selector, l, op, site1 => site2, field, boundaries)
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
