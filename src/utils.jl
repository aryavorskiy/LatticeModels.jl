using SparseArrays, StaticArrays
@static if VERSION < v"1.8"
    allequal(seq) = all(s == first(seq) for s in seq)
end

struct SparseMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    SparseMatrixBuilder{T}(sz) where T = new{T}(sz, [], [], [])
    SparseMatrixBuilder{T}(sz...) where T = SparseMatrixBuilder{T}(sz)
end
Base.size(smb::SparseMatrixBuilder) = smb.size
materialize(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)

Base.@propagate_inbounds function increment!(matrix::Matrix, rhs, idx1, idx2)
    matrix[idx1, idx2] .= rhs .+ @view matrix[idx1, idx2]
    return
end
Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs, idx1, idx2)
    for i in eachindex(idx1), j in eachindex(idx2)
        v = rhs[i, j]
        iszero(v) && continue
        push!(builder.Is, idx1[i])
        push!(builder.Js, idx2[j])
        push!(builder.Vs, v)
    end
end
function one_hot end
one_hot(index, ::Val{N}) where N = Int.(SVector{N}(1:N) .== index)
one_hot(index, N) = one_hot(index, Val(N))

const LenType{N} = Union{SVector{N}, NTuple{N, Any}}
@inline @generated function dot_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    Expr(:call, :+, [:(p1[$i] * p2[$i]) for i in 1:min(M, N)]...)
end
@inline @generated function add_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    mi = min(M, N)
    Expr(:call, :(SVector{M}), [:(p1[$i] + p2[$i]) for i in 1:mi]...,
        [M > N ? :(p1[$i]) : :(p2[$i]) for i in mi+1:M]...)
end
@inline @generated function mm_assuming_zeros(m, v::LenType{N}) where N
    Expr(:call, :+, [:(m[:, $i] * v[$i]) for i in 1:N]...)
end

function dims end
function lattice end
