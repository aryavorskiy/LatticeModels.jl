using SparseArrays, StaticArrays
@static if VERSION < v"1.8"
    allequal(seq) = all(s == first(seq) for s in seq)
end

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
