using LinearAlgebra, Statistics, Logging
import Base: length, getindex, show, copy, ==, zero

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

function show(io::IO, m::MIME"text/plain", b::Basis)
    println(io, "Basis with $(b.internal_dim)-dimensional internal phase space")
    print(io, "on ")
    show(io, m, b.lattice)
end

@doc """
    LatticeArray{LT, MT} where {LT<:Lattice, MT<:AbstractArray}

A wrapper object for array representing a wave function or linear operator.
Stores information about its basis to perform lattice checks.
"""
struct LatticeArray{LT<:Lattice,MT<:AbstractArray}
    basis::Basis{LT}
    operator::MT
    function LatticeArray(basis::Basis{LT}, operator::MT) where {LT<:Lattice,MT<:AbstractArray}
        !all((ax in (1, length(basis))) for ax in size(operator)) && error("inconsistent vector/matrix size")
        new{LT,MT}(basis, operator)
    end
end

const LatticeVector{LT,T} = LatticeArray{LT,T} where {LT<:Lattice,T<:AbstractVector}
const LatticeOperator{LT,T} = LatticeArray{LT,T} where {LT<:Lattice,T<:AbstractMatrix}

"""
    LatticeOperator{LT, MT}

The same as `LatticeArray{LT, MT}` where `MT<:AbstractMatrix`.

---
    LatticeOperator(uniform_scaling, basis)

Creates a `LatticeOperator` representation of a `UniformScaling` operator on given basis.
"""
function LatticeOperator(op::UniformScaling, bas::Basis)
    N = dims_internal(bas)
    m = Matrix(op, N, N)
    diag_operator(lattice(bas), m)
end

size(la::LatticeArray) = size(la.operator)
lattice(la::LatticeArray) = lattice(la.basis)
basis(la::LatticeArray) = la.basis
dims_internal(x) = dims_internal(basis(x))

@inline _ranges(is::Tuple, N::Int) = _ranges((), is, N)
@inline _ranges(rngs::Tuple, ::Tuple{}, N::Int) = rngs
@inline _ranges(rngs::Tuple, is::Tuple, N::Int) = _ranges(rngs, is[1], Base.tail(is), N)
@inline _ranges(rngs::Tuple, i::Int, is::Tuple, N::Int) =
    _ranges((rngs..., N*(i-1)+1:N*i), is, N)
@inline _ranges(rngs::Tuple, ::Colon, is::Tuple, N::Int) =
    _ranges((rngs..., :), is, N)
getindex(lo::LatticeArray, is::Int...) = lo.operator[_ranges(is, dims_internal(lo))...]
setindex!(lo::LatticeArray, val, is::Int...) =
    (lo.operator[_ranges(is, dims_internal(lo))...] = val)

==(lvm1::LatticeArray, lvm2::LatticeArray) = (lvm1.basis == lvm2.basis) && (lvm1.operator == lvm2.operator)

_typename(::LatticeVector) = "LatticeVector"
_typename(::LatticeOperator) = "LatticeOperator"
_typename(::LatticeArray) = "LatticeArray"
function show(io::IO, m::MIME"text/plain", la::LatticeArray{LT,MT} where {LT}) where {MT}
    println(io, join(size(la), "×"), " ", _typename(la), " with inner type $MT")
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
_zero_on_basis(bas::Basis) = _zero_on_basis(lattice(bas), dims_internal(bas))
function _zero_on_basis(l::Lattice, tp::TensorProduct)
    l != lattice(tp) &&
        error("lattice mismatch:\n$l\n$(lattice(tp))")
    zero(tp)
end

@inline _get_matrix_value(f::Function, l::Lattice, site::LatticeSite, ::Int) = f(site, site_coords(l, site))
@inline _get_matrix_value(m::AbstractMatrix, ::Lattice, ::LatticeSite, ::Int) = m
@inline _get_matrix_value(tp::TensorProduct, ::Lattice, ::LatticeSite, i::Int) = tp.lattice_value.values[i] * tp.matrix
function _diag_operator!(lop::LatticeOperator, op_object)
    l = lattice(lop)
    i = 1
    try
        for site in l
            lop[i, i] += _get_matrix_value(op_object, l, site, i)
            i += 1
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
    diag_operator(site_fun, lattice)

Creates a diagonal operator by applying the `site_fun` function to each site of `lattice`.
The `site_fun` function must accept a `LatticeSite` and its coordinates and return a matrix
which will be an operator affecting the internal state of the site.
"""
diag_operator(f::Function, l::Lattice) = _diag_operator!(_zero_on_basis(l, f), f)

"""
    diag_operator(lattice, matrix)

Creates a diagonal operator which affects only the internal state the same way on every site.
`matrix` is an `AbstractMatrix` representing the linear operator on the internal space.

Note that the matrix of the output `LatticeOperator` will be similar to `matrix`:
for instance, if `matrix` is sparse, so will be the output.
"""
diag_operator(l::Lattice, m::AbstractMatrix) = _diag_operator!(_zero_on_basis(l, m), m)

"""
    diag_operator(site_fun, basis)

Creates a diagonal operator which affects only the lattice space.
The `site_fun` function must accept a `LatticeSite` and its coordinates and return a number
which will be the diagonal element of the operator in lattice space.
"""
function diag_operator(lf::Function, bas::Basis)
    N = dims_internal(bas)
    eye = Matrix(I, N, N)
    _diag_operator!(_zero_on_basis(bas), (site, crd) -> lf(site, crd)::Number * eye)
end

"""
    diag_operator(basis, lattice_value)

Creates a diagonal operator which affects only the lattice space.
The `lattice_value` argument must be a `LatticeValue` storing diagonal elements of the operator in lattice space.
"""
function diag_operator(bas::Basis, lv::LatticeValue{<:Number})
    N = dims_internal(bas)
    eye = Matrix(I, N, N)
    _diag_operator!(_zero_on_basis(bas), TensorProduct(lv, eye))
end

"""
    coord_operators(basis)

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
    return xyz_operators
end

"""
    coord_operators(lattice, ndims)

The same as `coord_operators(Basis(lattice, ndims))`.
"""
coord_operators(l::Lattice, N::Int) = coord_operators(Basis(l, N))

"""
    diag_aggregate(function, lattice_operator)

Creates a `LatticeValue` where a site maps to the result of `function` on the matrix
of the operator narrowed to that site.
"""
diag_aggregate(f::Function, lo::LatticeArray) =
    LatticeValue(lattice(lo), [f(lo[i, i]) for i in 1:length(lattice(lo))])

"""
    ptrace(lattice_operator)

Same as `diag_aggregate(tr, lattice_operator)`
"""
ptrace(lo::LatticeArray) = diag_aggregate(tr, lo)

"""
    site_density(lattice_vector)
    site_density(lattice_operator)

A convenience function to find local density for wave functions (represented by `lattice_vector`)
and density matrices(represented by `lattice_operator`).

For wave functions it yields the total probability of the particle of being on every site.
For density matrices the partial trace is returned.
The return type is `LatticeValue` for both types of arguments.
"""
site_density(lo::LatticeOperator) = real.(ptrace(lo))
function site_density(lv::LatticeVector)
    l = lattice(lv)
    LatticeValue(l, [norm(lv[i])^2 for i in 1:length(l)])
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
    _unwrap_wlattice(f, el.basis, (checked_args..., el.operator), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)

@inline _unwrap_nolattice(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)
@inline function _unwrap_nolattice(f::Function, checked_args::Tuple, el::LatticeArray, args::Tuple; kw...)
    @warn "avoid using lattice operators and matrices in one function call"
    _unwrap_wlattice(f, el.basis, (checked_args..., el.operator), args; kw...)
end
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_nolattice(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap_nolattice(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, (checked_args..., el), args; kw...)
@inline function _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...)
    @warn "avoid using lattice operators and matrices in one function call"
    _unwrap_wlattice(f, basis, (checked_args..., el), args; kw...)
end
@inline function _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, el::LatticeArray, args::Tuple; kw...)
    el.basis != basis && throw(ArgumentError("basis mismatch:\n$(el.basis)\n$basis"))
    _unwrap_wlattice(f, basis, (checked_args..., el.operator), args; kw...)
end
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, checked_args, args[begin], Base.tail(args); kw...)
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, ::Tuple{}; kw...) =
    _make_wrapper(f(checked_args...; kw...), basis)

LatticeSummable = Union{LatticeArray,UniformScaling}
import Base: +, -, *, /, ^, adjoint, copy
@inline +(los::LatticeSummable...) = _unwrap(+, los)
@inline -(lo::LatticeArray) = _unwrap(-, (lo,))
@inline -(lo1::LatticeSummable, lo2::LatticeSummable) = _unwrap(-, (lo1, lo2))
@inline *(los::LatticeArray...) = _unwrap(*, los)
@inline *(num::Number, lo::LatticeArray) = _unwrap(*, (num, lo))
@inline *(lo::LatticeArray, num::Number) = _unwrap(*, (lo, num))
@inline /(lo::LatticeArray, num::Number) = _unwrap(/, (lo, num))
@inline ^(lo::LatticeArray, num::Number) = _unwrap(^, (lo, num))
@inline adjoint(lo::LatticeArray) = _unwrap(adjoint, (lo,))
@inline copy(lo::LatticeArray) = _unwrap(copy, (lo,))

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
    @on_lattice

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
xexpy == @on_lattice X * exp(Y) + Y     # true
```
"""
macro on_lattice(expr)
    we = _wrap_smart!(expr)
    return we
end
