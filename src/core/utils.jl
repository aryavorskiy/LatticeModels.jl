using SparseArrays, StaticArrays

@generated function one_hot(indices, ::Val{N}) where N
    Expr(:call, :(SVector{N}), (:(Int($v in indices)) for v in 1:N)...)
end
one_hot(indices, N) = one_hot(indices, Val(N))

const LenType{N} = Union{SVector{N}, NTuple{N, Any}}
@inline @generated function dot_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    Expr(:call, :+, Expr[:(p1[$i] * p2[$i]) for i in 1:min(M, N)]...)
end
@inline @generated function add_assuming_zeros(p1::LenType{M}, p2::LenType{N}) where {M,N}
    out_expr = Expr(:call, :(SVector{M}))
    for i in 1:M
        if i ≤ N
            push!(out_expr.args, :(p1[$i] + p2[$i]))
        else
            push!(out_expr.args, :(p1[$i]))
        end
    end
    return out_expr
end

@generated cartesian_indices(depth::Int, ::Val{M}) where M = quote
    CartesianIndex($((:(-depth) for _ in 1:M)...)):CartesianIndex($((:depth for _ in 1:M)...))
end

# Size check
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
            "invalid size of `$var`; expected square matrix, got size ", join(size(arr), "×"))))
    end
end

macro check_size(arr::Symbol, siz)
    quote
        check_size($(esc(arr)), $(esc(siz)), $(QuoteNode(arr)))
    end
end

# Filter tuples at compile-time
skiptype(T::Type, args::Tuple) = skiptype(T, (), args)
skiptype(::Type{T}, pre_args::Tuple, post_args::Tuple) where T =
    first(post_args) isa T ?
    skiptype(T, pre_args, Base.tail(post_args)) :
    skiptype(T, (pre_args..., first(post_args)), Base.tail(post_args))
skiptype(::Type, pre_args::Tuple, ::Tuple{}) = pre_args

# Indentation
getindent(io::IO) = "  " ^ get(io, :indent, 0)
addindent(io::IO, n::Int, pairs...) = IOContext(io, :indent => get(io, :indent, 0) + n, pairs...)
addindent(io::IO, pairs...) = addindent(io, 1, pairs...)

# Compact display if requested or if showing a set
requires_compact(io::IO) = get(io, :compact, false) || get(io, :SHOWN_SET, nothing) !== nothing

# Format a number with a noun properly
function fmtnum(n::Int, noun::String, singular::String, plural::String)
    suffix = (n % 10 == 1 && n % 100 != 11) ? singular : plural
    return "$n $noun$suffix"
end
fmtnum(n::Int, noun::String) = fmtnum(n, noun, "", "s")
fmtnum(any, args...) = fmtnum(length(any), args...)
