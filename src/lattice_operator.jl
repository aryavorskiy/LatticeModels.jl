using LinearAlgebra, Statistics, Logging
import Base: length, getindex, show, copy, ==, zero

struct Basis{LT<:Lattice}
    lattice::LT
    internal_dim::Int
end

length(b::Basis) = length(b.lattice) * b.internal_dim
dims_internal(b::Basis) = b.internal_dim

function show(io::IO, m::MIME"text/plain", b::Basis)
    println(io, "Basis with $(b.internal_dim)-dimensional internal phase space")
    print(io, "on ")
    show(io, m, b.lattice)
end

struct LatticeArray{LT<:Lattice,MT<:AbstractArray}
    basis::Basis{LT}
    operator::MT
    function LatticeArray(basis::Basis{LT}, operator::MT) where {LT<:Lattice,MT<:AbstractArray}
        !all(size(operator) .== length(basis)) && error("inconsistent vector/matrix size")
        new{LT,MT}(basis, operator)
    end
end

const LatticeVector{LT,T} = LatticeArray{LT,T} where {LT<:Lattice,T<:AbstractVector}
const LatticeOperator{LT,T} = LatticeArray{LT,T} where {LT<:Lattice,T<:AbstractMatrix}

function LatticeOperator(op::UniformScaling, bas::Basis)
    N = bas.internal_dim
    m = Matrix(op, N, N)
    diag_operator(bas.lattice, m)
end

size(lv::LatticeArray) = size(lv.operator)
dims_internal(lv::LatticeArray) = lv.basis.internal_dim

@inline _ranges(is::Tuple, N::Int) = _ranges((), is, N)
@inline _ranges(rngs::Tuple, ::Tuple{}, N::Int) = rngs
@inline _ranges(rngs::Tuple, is::Tuple, N::Int) = _ranges(rngs, is[1], Base.tail(is), N)
@inline _ranges(rngs::Tuple, i::Int, is::Tuple, N::Int) =
    _ranges((rngs..., N*(i-1)+1:N*i), is, N)
getindex(lo::LatticeArray, is::Int...) = lo.operator[_ranges(is, dims_internal(lo))...]
setindex!(lo::LatticeArray, val, is::Int...) =
    (lo.operator[_ranges(is, dims_internal(lo))...] = val)

==(lvm1::LatticeArray, lvm2::LatticeArray) = (lvm1.basis == lvm2.basis) && (lvm1.operator == lvm2.operator)

function show(io::IO, m::MIME"text/plain", lv::LatticeArray{LT,MT}) where {LT,MT<:AbstractMatrix}
    println(io, join(size(lv), "×") * " LatticeMatrix with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeArray{MT}) where {MT<:AbstractVector}
    println(io, "$(length(lv.operator))-length LatticeVector with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

struct TensorProduct{LVT<:LatticeValue{<:Number},MT<:AbstractMatrix}
    lattice_value::LVT
    matrix::MT
end

dims_internal(tp::TensorProduct) = size(tp.matrix)[1]
basis(tp::TensorProduct) = Basis(tp.lattice_value.lattice, dims_internal(tp))
zero(tp::TensorProduct) = _zero_on_basis(tp.lattice_value.lattice, tp.matrix)
copy(tp::TensorProduct) = materialize(tp)
⊗(lv::LatticeValue, m::Matrix) = copy(TensorProduct(lv, m))
⊗(m::Matrix, lv::LatticeValue) = copy(TensorProduct(lv, m))

function _zero_on_basis(l::Lattice, m::AbstractMatrix)
    N = size(m)[1]
    LatticeArray(Basis(l, N),
        zero(similar(m, ComplexF64, (N * length(l), N * length(l)))))
end
function _zero_on_basis(l::Lattice, lf::Function)
    _zero_on_basis(l, lf(first(l), coords(l, first(l))))
end
_zero_on_basis(l::Lattice, N::Int) = LatticeArray(Basis(l, N),
    zeros(ComplexF64, N * length(l), N * length(l)))
_zero_on_basis(bas::Basis) = _zero_on_basis(bas.lattice, bas.internal_dim)
function _zero_on_basis(l::Lattice, tp::TensorProduct)
    l != tp.lattice_value.lattice &&
        error("lattice mismatch:\n$l\n$(tp.lattice_value.lattice)")
    zero(tp)
end

@inline _get_matrix_value(f::Function, l::Lattice, site::LatticeSite, ::Int) = f(site, coords(l, site))
@inline _get_matrix_value(m::AbstractMatrix, ::Lattice, ::LatticeSite, ::Int) = m
@inline _get_matrix_value(tp::TensorProduct, ::Lattice, ::LatticeSite, i::Int) = tp.lattice_value.values[i] * tp.matrix
function _diag_operator!(lop::LatticeOperator, op_object)
    l = lop.basis.lattice
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
diag_operator(f::Function, l::Lattice) = _diag_operator!(_zero_on_basis(l, f), f)
diag_operator(l::Lattice, m::AbstractMatrix) = _diag_operator!(_zero_on_basis(l, m), m)
function diag_operator(lf::Function, bas::Basis)
    N = bas.internal_dim
    eye = Matrix(I, N, N)
    _diag_operator!(_zero_on_basis(bas), (site, crd) -> lf(site, crd)::Number * eye)
end
function diag_operator(bas::Basis, lv::LatticeValue{<:Number})
    N = bas.internal_dim
    eye = Matrix(I, N, N)
    _diag_operator!(_zero_on_basis(bas), TensorProduct(lv, eye))
end

function coord_operators(bas::Basis)
    N = bas.internal_dim
    d = dims(bas.lattice)
    i = 1
    eye = Matrix(I, N, N)
    xyz_operators = [LatticeArray(bas, op_mat) for op_mat in
                     eachslice(zeros(length(bas), length(bas), d), dims=3)]
    for site in bas.lattice
        crd = coords(bas.lattice, site)
        for j in 1:d
            xyz_operators[j][i, i] = crd[j] * eye
        end
        i += 1
    end
    return xyz_operators
end

coord_operators(l::Lattice, N::Int) = coord_operators(Basis(l, N))

diag_aggregate(f::Function, lo::LatticeArray) =
    LatticeValue(lo.basis.lattice, [f(lo[i, i]) for i in 1:length(lo.basis.lattice)])

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

macro on_lattice(expr)
    we = _wrap_smart!(expr)
    return we
end
