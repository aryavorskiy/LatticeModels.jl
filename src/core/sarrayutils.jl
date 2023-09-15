using SparseArrays, StaticArrays

@generated function one_hot(indices, ::Val{N}) where N
    Expr(:call, :(SVector{N}), (:(Int($v in indices)) for v in 1:N)...)
end
one_hot(indices, N) = one_hot(indices, Val(N))

const LenType{N} = Union{SVector{N}, NTuple{N, Any}}
@inline @generated function dot_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    Expr(:call, :+, [:(p1[$i] * p2[$i]) for i in 1:min(M, N)]...)
end
@inline @generated function add_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    mi = min(M, N)
    Expr(:call, :(SVector{M}), [:(p1[$i] + p2[$i]) for i in 1:mi]...,
        [M > N ? :(p1[$i]) : :(p2[$i]) for i in mi+1:M]...)
end
@inline @generated function mm_assuming_zeros(m::SMatrix{M}, v::LenType{N}) where {M, N}
    if N == 0
        return :(zero(SVector{$M}))
    else
        Expr(:call, :+, [:(m[:, $i] * v[$i]) for i in 1:N]...)
    end
end

const Nullable{T} = Union{Nothing,T}
