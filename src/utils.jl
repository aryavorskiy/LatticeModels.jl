using SparseArrays, StaticArrays
import Base: size, materialize
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
size(smb::SparseMatrixBuilder) = smb.size
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
