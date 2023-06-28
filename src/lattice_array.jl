using LinearAlgebra, Logging
import Base: ==

abstract type Basis end
basis(b::Basis) = b

"""
    check_basis_match(b1, b2)

Checks if the two objects are defined on the same basis, throws an error if not.
`b1` and `b2` can be of any type that defines a method for the `basis` function.
"""
function check_basis_match(b1, b2)
    basis(b1) != basis(b2) &&
        throw(ArgumentError("basis mismatch:\n$(repr("text/plain", basis(b1)))\n$(repr("text/plain", basis(b2)))"))
end
function to_slice end

@doc """
    LatticeArray{AT, BT, N}

A wrapper object for array representing a wavefunction or linear operator.
Stores information about its basis to perform lattice checks.
"""
struct LatticeArray{AT,BT,N}
    basis::BT
    array::AT
    function LatticeArray{N}(basis::BT, array::AT) where {N, BT<:Basis, AT}
        !all((ax in (1, length(basis))) for ax in size(array)) &&
            throw(DimensionMismatch("array has size $(size(array)), basis has length $(length(basis))"))
        new{AT,BT,N}(basis, array)
    end
    LatticeArray(basis, arr::AbstractArray{_T, N} where _T) where N =
        LatticeArray{N}(basis, arr)
    LatticeArray(basis, smb::SparseMatrixBuilder) =
        LatticeArray{2}(basis, smb)
end

const LatticeVector{VT,BT} = LatticeArray{VT,BT,1}
const LatticeOperator{MT,BT} = LatticeArray{MT,BT,2}

"""
    LatticeOperator{MT, BT}

The same as `LatticeArray{MT, BT, 2}` where `MT<:AbstractMatrix`.

---
    LatticeOperator(uniform_scaling::UniformScaling, basis::Basis)

Creates a `LatticeOperator` representation of a `UniformScaling` operator on given basis.
For example, `LatticeOperator(LinearAlgebra.I, basis)` yields an identity operator on `basis` basis.
"""
LatticeOperator(bas::Basis, op::UniformScaling) =
    LatticeArray(bas, Matrix(op, length(bas), length(bas)))

basis(la::LatticeArray) = la.basis
dims_internal(x) = dims_internal(basis(x))
lattice(x) = lattice(basis(x))

Base.zero(la::LatticeArray) = LatticeArray(basis(la), zero(la))
Base.copy(la::LatticeArray) = LatticeArray(basis(la), copy(la.array))
Base.size(la::LatticeArray) = size(la.array)
==(lvm1::LatticeArray, lvm2::LatticeArray) = (lvm1.basis == lvm2.basis) && (lvm1.array == lvm2.array)

@inline _to_indices(is::Tuple, b::Basis) = _to_indices((), is, b)
@inline _to_indices(rngs::Tuple, ::Tuple{}, ::Basis) = rngs
@inline _to_indices(rngs::Tuple, is::Tuple, b::Basis) = _to_indices(rngs, is[1], Base.tail(is), b)
@inline _to_indices(rngs::Tuple, i, is::Tuple, b::Basis) =
    _to_indices((rngs..., to_slice(b, i)), is, b)
Base.getindex(la::LatticeArray, is::Vararg{Any}) = la.array[_to_indices(is, basis(la))...]
Base.view(la::LatticeArray, is::Vararg{Any}) = view(la.array, _to_indices(is, basis(la))...)
Base.setindex!(la::LatticeArray, val, is::Vararg{Any}) =
    (la.array[_to_indices(is, basis(la))...] = val)
increment!(la::LatticeArray, rhs, is::Vararg{Any}) =
    increment!(la.array, rhs, _to_indices(is, basis(la))...)

_typename(::LatticeVector) = "LatticeVector"
_typename(::LatticeOperator) = "LatticeOperator"
_typename(::LatticeArray) = "LatticeArray"
function Base.show(io::IO, m::MIME"text/plain", la::LatticeArray{AT}) where {AT}
    print(io, join(size(la), "×"))
    AT<:LatticeVector && print(io, "-element")
    println(io, " ", _typename(la), " with inner type $AT")
    print(io, "on ")
    show(io, m, la.basis)
end

"""
    diag_reduce(f, lattice_operator::LatticeOperator)

Creates a `LatticeValue` where a site maps to the result of `f` on the matrix
of the operator narrowed to that site.
"""
diag_reduce(f::Function, lo::LatticeOperator) =
    LatticeValue(lattice(lo), [f(lo[i, i]) for i in 1:length(lattice(lo))])

"""
    ptrace(lattice_operator::LatticeOperator, space)

Calculates a matrix for the partial trace of given operator.
`space` argument must take one of two values:
- `:lattice` for taking the partial trace over the lattice space.
- `:internal` for the same over the internal space.
"""
function ptrace(lo::LatticeOperator, space::Symbol)
    N = dims_internal(lo)
    if space === :lattice
        sum(@views lo[i, i] for i in 1:length(lattice(lo)))
    elseif space === :internal
        blen = length(basis(lo))
        sum(@views lo.array[i:N:blen, i:N:blen] for i in 1:N)
    else
        throw(ArgumentError("unsupported value '$space' of 'space' argument"))
    end
end

"""
    site_density(lattice_vector::LatticeVector)
    site_density(lattice_operator::LatticeOperator)

A convenience function to find local density for wavefunctions (represented by `lattice_vector`)
and density matrices (represented by `lattice_operator`).

Returns a LatticeValue representing the total probability of the particle of being on every site.
"""
site_density(lo::LatticeOperator) =
    LatticeValue(lattice(lo), vec(real.(sum(reshape(diag(lo.array), (dims_internal(lo), :)), dims=1))))
function site_density(lv::LatticeVector)
    l = lattice(lv)
    LatticeValue(l, [sum(abs2, lv[i]) for i in 1:length(l)])
end

@inline _make_wrapper(op, ::Nothing) = op
@inline _make_wrapper(op::Number, ::Basis) = op
@inline _make_wrapper(op::Any, b::Basis) = LatticeArray(b, op)

@inline _unwrap_from_macro(f, args...; kw...) = _unwrap(f, args; kw...)
@inline _unwrap(T::Type, args::Tuple; kw...) = T(args...; kw...)

@inline _unwrap(f::Function, args::Tuple; kw...) = _unwrap(f, (), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap(f, (checked_args..., el), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap(f::Function, checked_args::Tuple, el::LatticeArray, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis(el), (checked_args..., el.array), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)

@inline _unwrap_nolattice(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)
@inline function _unwrap_nolattice(f::Function, checked_args::Tuple, el::LatticeArray, args::Tuple; kw...)
    @warn "avoid using lattice arrays and unwrapped arrays in one function call"
    _unwrap_wlattice(f, basis(el), (checked_args..., el.array), args; kw...)
end
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_nolattice(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, (checked_args..., el), args; kw...)
@inline function _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...)
    @warn "avoid using lattice arrays and unwrapped arrays in one function call"
    _unwrap_wlattice(f, basis, (checked_args..., el), args; kw...)
end
@inline function _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, el::LatticeArray, args::Tuple; kw...)
    check_basis_match(el, basis)
    _unwrap_wlattice(f, basis, (checked_args..., el.array), args; kw...)
end
@inline _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, checked_args, args[begin], Base.tail(args); kw...)
@inline _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, ::Tuple{}; kw...) =
    _make_wrapper(f(checked_args...; kw...), basis)

LatticeSummable = Union{LatticeArray,UniformScaling}
import Base: +, -, *, /, \, ^, adjoint, exp, inv
@inline +(ls::LatticeSummable, lss::LatticeSummable...) = _unwrap(+, (ls, lss...))
@inline -(lo1::LatticeSummable, lo2::LatticeSummable) = _unwrap(-, (lo1, lo2))
@inline *(la::LatticeArray, las::LatticeArray...) = _unwrap(*, (la, las...))
for f in (:*, :/, :\, :^)
    @eval @inline ($f)(la::LatticeArray, num::Number) = _unwrap(($f), (la, num))
    @eval @inline ($f)(num::Number, la::LatticeArray) = _unwrap(($f), (num, la))
end
for f in (:adjoint, :exp, :inv, :-)
    @eval @inline ($f)(la::LatticeArray) = _unwrap(($f), (la,))
end

import LinearAlgebra: dot
@inline dot(lv1::LatticeVector, lv2::LatticeVector) = _unwrap(dot, (lv1, lv2))
@inline dot(lv1::LatticeVector, A::LatticeOperator, lv2::LatticeVector) = _unwrap(dot, (lv1, A, lv2))

_wrap_smart!(expr::Any) = expr
function _wrap_smart!(expr::Expr)
    Meta.isexpr(expr, :escape) && error("do not use @on_lattice macro in other macros")
    local _begin = 1
    if Meta.isexpr(expr, :call)
        insert!(expr.args, 1, :(_unwrap_from_macro))
        _begin = 2
    elseif Meta.isexpr(expr, (:function, :->, :kw))
        _begin = 2
    end
    for i in _begin:length(expr.args)
        if expr.args[i] isa Symbol
            expr.args[i] = :($(esc(expr.args[i])))
        else
            _wrap_smart!(expr.args[i])
        end
    end
    return expr
end

"""
    @on_lattice expression

Replaces all `LatticeArray`s in subsequent function calls with actual arrays stored inside them.
Throws an error if lattice arrays in one function call are defined on different lattices,
shows a warning if a lattice array is used in one call with a normal array.

## Example
```julia
l = SquareLattice(10, 10)
bas = Basis(l, 2)
X, Y = coord_operators(bas)
xexpypy = diag_operator(bas) do (x, y)
    x * exp(y) + y
end
xexpypy == @on_lattice X * exp(Y) + Y     # true
```
"""
macro on_lattice(expr)
    _wrap_smart!(expr)
end
