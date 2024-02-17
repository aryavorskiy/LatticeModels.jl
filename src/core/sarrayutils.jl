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

@generated cartesian_indices(depth::Int, ::Val{M}) where M = quote
    CartesianIndex($((:(-depth) for _ in 1:M)...)):CartesianIndex($((:depth for _ in 1:M)...))
end

function check_size(arr, expected_size::Tuple, var::Symbol)
    if size(arr) != expected_size
        throw(ArgumentError(string(
            "invalid ",
            arr isa AbstractVector ? "length" : "size",
            " of `$var`; expected ",
            join(expected_size, "×"),
            ", got ", join(size(arr), "×"))))
    end
end
check_size(arr, len::Int, var) = check_size(arr, (len,), var)
function check_size(arr, expected_size::Symbol, var::Symbol)
    if expected_size === :square
        length(size(arr)) == 2 && size(arr, 1) == size(arr, 2) && return
        throw(ArgumentError(string(
            "invalid size of `$var`; expected square matrix, got ", join(size(arr), "×"))))
    end
end

macro check_size(arr::Symbol, siz)
    quote
        check_size($(esc(arr)), $(esc(siz)), $(QuoteNode(arr)))
    end
end
