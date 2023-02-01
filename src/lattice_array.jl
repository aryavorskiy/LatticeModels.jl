using LinearAlgebra, Statistics, Logging
import Base: length, getindex, view, show, copy, ==, zero

"""
    Basis{LT} where {LT<:Lattice}

A basis on a lattice with some number of internal states on each site.
Fields:
- `lattice`: the [`Lattice`](@ref) of the basis
- `internal_dim`: the number of internal states on each site
"""
struct Basis{LT<:Lattice}
    lattice::LT
    internal_dim::Int
end
lattice(b::Basis) = b.lattice
dims_internal(b::Basis) = b.internal_dim
length(b::Basis) = length(lattice(b)) * dims_internal(b)
==(b1::Basis, b2::Basis) = b1.internal_dim == b2.internal_dim && b1.lattice == b2.lattice

basis(b::Basis) = b
function check_basis_match(b1, b2)
    basis(b1) != basis(b2) &&
        throw(ArgumentError("basis mismatch:\n$(_s(basis(b1)))\n$(_s(basis(b2)))"))
end

function show(io::IO, m::MIME"text/plain", b::Basis)
    println(io, "Basis with $(b.internal_dim)-dimensional internal phase space")
    print(io, "on ")
    show(io, m, b.lattice)
end

@doc """
    LatticeArray{AT, LT, N}

A wrapper object for array representing a wave function or linear operator.
Stores information about its basis to perform lattice checks.
"""
struct LatticeArray{AT,LT,N}
    basis::Basis{LT}
    array::AT
    function LatticeArray(basis::Basis{LT}, array::AT) where {LT<:Lattice,AT<:AbstractArray{_T,N} where _T} where {N}
        !all((ax in (1, length(basis))) for ax in size(array)) &&
            throw(DimensionMismatch("array has size $(size(array)), basis has length $(length(basis))"))
        new{AT,LT,N}(basis, array)
    end
end

const LatticeVector{VT,LT} = LatticeArray{VT,LT,1}
const LatticeOperator{MT,LT} = LatticeArray{MT,LT,2}

"""
    LatticeOperator{MT, LT}

The same as `LatticeArray{MT, LT, 2}` where `MT<:AbstractMatrix`.

---
    LatticeOperator(uniform_scaling::UniformScaling, basis::Basis)

Creates a `LatticeOperator` representation of a `UniformScaling` operator on given basis.
For example, `LatticeOperator(LinearAlgebra.I, basis)` yields an identity operator on `basis` basis.
"""
function LatticeOperator(bas::Basis, op::UniformScaling)
    N = dims_internal(bas)
    m = Matrix(op, N, N)
    diag_operator(lattice(bas), m)
end

size(la::LatticeArray) = size(la.array)
basis(la::LatticeArray) = la.basis
dims_internal(x) = dims_internal(basis(x))
lattice(x) = lattice(basis(x))

@inline _ranges(is::Tuple, l::Lattice, N::Int) = _ranges((), is, l, N)
@inline _ranges(rngs::Tuple, ::Tuple{}, ::Lattice, N::Int) = rngs
@inline _ranges(rngs::Tuple, is::Tuple, l::Lattice, N::Int) = _ranges(rngs, is[1], Base.tail(is), l, N)
@inline _ranges(rngs::Tuple, i::Int, is::Tuple, l::Lattice, N::Int) =
    _ranges((rngs..., N*(i-1)+1:N*i), is, l, N)
@inline _ranges(rngs::Tuple, site::LatticeSite, is::Tuple, l::Lattice, N::Int) =
    _ranges((rngs..., site_index(l, site)), is, l, N)
@inline _ranges(rngs::Tuple, ::Colon, is::Tuple, l::Lattice, N::Int) =
    _ranges((rngs..., :), is, l, N)
getindex(la::LatticeArray, is::Vararg{Any}) = la.array[_ranges(is, lattice(la), dims_internal(la))...]
Base.view(la::LatticeArray, is::Vararg{Any}) = view(la.array, _ranges(is, lattice(la), dims_internal(la))...)
setindex!(la::LatticeArray, val, is::Vararg{Any}) =
    (la.array[_ranges(is, lattice(la), dims_internal(la))...] = val)

==(lvm1::LatticeArray, lvm2::LatticeArray) = (lvm1.basis == lvm2.basis) && (lvm1.array == lvm2.array)

_typename(::LatticeVector) = "LatticeVector"
_typename(::LatticeOperator) = "LatticeOperator"
_typename(::LatticeArray) = "LatticeArray"
function show(io::IO, m::MIME"text/plain", la::LatticeArray{AT}) where {AT}
    print(io, join(size(la), "×"))
    AT<:LatticeVector && print(io, "-element")
    println(io, " ", _typename(la), " with inner type $AT")
    print(io, "on ")
    show(io, m, la.basis)
end

"""
    TensorProduct{LVT, MT} where {LVT<:LatticeValue{<:Number}, MT<:AbstractMatrix}

A lazy representation of an operator as a tensor product of two distinct phase spaces.
One affects only the internal space, the other - only the lattice space.

The `lattice_value ⊗ matrix` notation computes the value of the `TensorProduct` eagerly,
which means that the result will be a `LatticeOperator`.
However, in the `@hamiltonian` macro lazy computation is forced.
"""
struct TensorProduct{LVT<:LatticeValue{<:Number},MT<:AbstractMatrix}
    lattice_value::LVT
    matrix::MT
end

dims_internal(tp::TensorProduct) = size(tp.matrix)[1]
lattice(tp::TensorProduct) = lattice(tp.lattice_value)
basis(tp::TensorProduct) = Basis(lattice(tp), dims_internal(tp))
zero(tp::TensorProduct) = _zero_on_basis(lattice(tp), tp.matrix)
copy(tp::TensorProduct) = materialize(tp)
⊗(lv::LatticeValue, m::Matrix) = copy(TensorProduct(lv, m))
⊗(m::Matrix, lv::LatticeValue) = copy(TensorProduct(lv, m))

function _zero_on_basis(l::Lattice, m::AbstractMatrix)
    N = size(m)[1]
    LatticeArray(Basis(l, N),
        zero(similar(m, ComplexF64, (N * length(l), N * length(l)))))
end
function _zero_on_basis(l::Lattice, lf::Function)
    _zero_on_basis(l, lf(first(l), site_coords(l, first(l))))
end
_zero_on_basis(l::Lattice, N::Int) = LatticeArray(Basis(l, N),
    zeros(ComplexF64, N * length(l), N * length(l)))
function _zero_on_basis(l::Lattice, N::Int, MT::Type)
    LatticeArray(Basis(l, N),
        zero(similar(MT, (N * length(l), N * length(l)))))
end
_zero_on_basis(l::Lattice, N::Int, ::Type{Matrix{ComplexF64}}) = _zero_on_basis(l, N)
_zero_on_basis(bas::Basis) = _zero_on_basis(lattice(bas), dims_internal(bas))
function _zero_on_basis(l::Lattice, tp::TensorProduct)
    check_lattice_match(l, tp)
    zero(tp)
end

@inline _get_matrix_value(f::Function, l::Lattice, site::LatticeSite, ::Int, ::Matrix) = f(site, site_coords(l, site))
@inline _get_matrix_value(m::AbstractMatrix, ::Lattice, ::LatticeSite, ::Int, ::Matrix) = m
@inline _get_matrix_value(tp::TensorProduct, ::Lattice, ::LatticeSite, i::Int, ::Matrix) = tp.lattice_value.values[i] * tp.matrix
@inline _get_matrix_value(lv::LatticeValue, ::Lattice, ::LatticeSite, i::Int, eye::Matrix) = lv.values[i] * eye
@inline _get_matrix_value(n::Number, ::Lattice, ::LatticeSite, i::Int, eye::Matrix) = n * eye
function _diag_operator!(lop::LatticeOperator, op_object)
    N = dims_internal(lop)
    eye = Matrix(I, N, N)
    l = lattice(lop)
    try
        for (i, site) in enumerate(l)
            lop[i, i] = _get_matrix_value(op_object, l, site, i, eye) + @view lop[i, i]
        end
    catch e
        if e isa DimensionMismatch
            error("dimension mismatch")
        else
            rethrow()
        end
    end
    lop
end

materialize(tp::TensorProduct) = _diag_operator!(_zero_on_basis(basis(tp)), tp)

"""
    diag_operator(f, l::Lattice)

Creates a diagonal operator by applying the `f` function to each site of `l`.
`f` must accept a `LatticeSite` and its coordinate vector and return a matrix
which will be an operator affecting the internal state of the site.
"""
diag_operator(f::Function, l::Lattice) = _diag_operator!(_zero_on_basis(l, f), f)

"""
    diag_operator(f, bas::Basis)

Creates a diagonal operator by applying the `f` function to each site of the lattice of `bas`.
`f` must accept a `LatticeSite` and its coordinate vector and return a matrix
which will be an operator affecting the internal state of the site.
"""
diag_operator(f::Function, bas::Basis) = _diag_operator!(_zero_on_basis(bas), f)

"""
    diag_operator(lattice::Lattice, matrix::AbstractMatrix)

Creates a diagonal operator which affects only the internal state the same way on every site.
`matrix` is an `AbstractMatrix` representing the linear operator on the internal space.

Note that the matrix of the output `LatticeOperator` will be similar to `matrix`:
for instance, if `matrix` is sparse, so will be the output.
"""
diag_operator(l::Lattice, m::AbstractMatrix) = _diag_operator!(_zero_on_basis(l, m), m)

"""
    diag_operator(lv::LatticeValue, N::Int=1)

Creates a diagonal operator which affects only the lattice space.
The `lv` argument must be a `LatticeValue` storing diagonal elements of the operator in lattice space.
`N` is the number of internal degrees of freedom on each site.
"""
function diag_operator(lv::LatticeValue{T}, N::Int=1) where T<:Number
    _diag_operator!(_zero_on_basis(lattice(lv), N, Array{T}), lv)
end

"""
    coord_operators(basis::Basis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
function coord_operators(bas::Basis)
    N = dims_internal(bas)
    d = dims(bas.lattice)
    i = 1
    eye = Matrix(I, N, N)
    xyz_operators = [LatticeArray(bas, op_mat) for op_mat in
                     eachslice(zeros(length(bas), length(bas), d), dims=3)]
    for site in bas.lattice
        crd = site_coords(bas.lattice, site)
        for j in 1:d
            xyz_operators[j][i, i] = crd[j] * eye
        end
        i += 1
    end
    xyz_operators
end

"""
    coord_operators(lattice::Lattice, ndims::Tnt)

The same as `coord_operators(Basis(lattice, ndims))`.
"""
coord_operators(l::Lattice, N::Int) = coord_operators(Basis(l, N))

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

A convenience function to find local density for wave functions (represented by `lattice_vector`)
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
    @warn "avoid using lattice operators and matrices in one function call"
    _unwrap_wlattice(f, basis(el), (checked_args..., el.array), args; kw...)
end
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_nolattice(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, (checked_args..., el), args; kw...)
@inline function _unwrap_wlattice(f::Function, basis::Basis, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...)
    @warn "avoid using lattice operators and matrices in one function call"
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
import Base: +, -, *, /, \, ^, adjoint, copy, exp, inv
@inline +(los::LatticeSummable...) = _unwrap(+, los)
@inline -(lo1::LatticeSummable, lo2::LatticeSummable) = _unwrap(-, (lo1, lo2))
@inline *(los::LatticeArray...) = _unwrap(*, los)
for f in (:*, :/, :\, :^)
    @eval @inline ($f)(la::LatticeArray, num::Number) = _unwrap(($f), (la, num))
    @eval @inline ($f)(num::Number, la::LatticeArray) = _unwrap(($f), (num, la))
end
for f in (:adjoint, :copy, :exp, :inv, :-)
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
Throws an error if lattice operators in one function call are defined on different lattices,
shows a warning if a lattice array is used in one call with a normal array.

## Example
```julia
l = SquareLattice(10, 10)
bas = Basis(l, 2)
X, Y = coord_operators(bas)
xexpypy = diag_operator(bas) do site, (x, y)
    x * exp(y) + y
end
xexpypy == @on_lattice X * exp(Y) + Y     # true
```
"""
macro on_lattice(expr)
    _wrap_smart!(expr)
end
