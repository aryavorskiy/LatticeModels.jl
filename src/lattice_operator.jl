using LinearAlgebra, Statistics, Logging
import Base: length, getindex, show, ==

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

struct LatticeVecOrMat{LT<:Lattice,MT<:AbstractVecOrMat}
    basis::Basis{LT}
    operator::MT
    function LatticeVecOrMat(basis::Basis{LT}, operator::MT) where {LT<:Lattice,MT<:AbstractVecOrMat}
        @assert all(size(operator) .== length(basis))
        new{LT,MT}(basis, operator)
    end
end

LatticeOperator{T} = LatticeVecOrMat{<:Lattice,T} where {T<:AbstractMatrix}

size(lv::LatticeVecOrMat) = size(lv.operator)
dims_internal(lv::LatticeVecOrMat) = lv.basis.internal_dim

@inline _ranges(is::Tuple, N::Int) = _ranges((), is, N)
@inline _ranges(rngs::Tuple, ::Tuple{}, N::Int) = rngs
@inline _ranges(rngs::Tuple, is::Tuple, N::Int) = _ranges(rngs, is[1], Base.tail(is), N)
@inline _ranges(rngs::Tuple, i::Int, is::Tuple, N::Int) =
    _ranges((rngs..., N*(i-1)+1:N*i), is, N)
getindex(lo::LatticeVecOrMat, is::Int...) = lo.operator[_ranges(is, dims_internal(lo))...]
setindex!(lo::LatticeVecOrMat, val, is::Int...) =
    (lo.operator[_ranges(is, dims_internal(lo))...] = val)

==(lvm1::LatticeVecOrMat, lvm2::LatticeVecOrMat) = (lvm1.basis == lvm2.basis) && (lvm1.operator == lvm2.operator)

function show(io::IO, m::MIME"text/plain", lv::LatticeVecOrMat{LT, MT}) where {LT, MT<:AbstractMatrix}
    println(io, join(size(lv), "×") * " LatticeMatrix with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeVecOrMat{MT}) where {MT<:AbstractVector}
    println(io, "$(length(lv.operator))-length LatticeVector with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

function _zero_on_basis(l::Lattice, m::AbstractMatrix)
    N = size(m)[1]
    LatticeVecOrMat(Basis(l, N),
        zero(similar(m, ComplexF64, (N * length(l), N * length(l)))))
end
function _zero_on_basis(l::Lattice, f::Function)
    lf = _propagate_lattice_args(f, l)
    _zero_on_basis(l, lf(l, first(l)))
end
_zero_on_basis(l::Lattice, N::Int) = LatticeVecOrMat(Basis(l, N),
    zeros(ComplexF64, N * length(l), N * length(l)))
_zero_on_basis(bas::Basis) = _zero_on_basis(bas.lattice, bas.internal_dim)

convert_inner_type(MT::Type, lo::LatticeVecOrMat) =
    LatticeVecOrMat(lo.basis, convert(MT, lo.operator))

# TODO: check if convert_inner_type() works

function _diag_operator!(lop::LatticeOperator, lf::Function)
    l = lop.basis.lattice
    i = 1
    try
        for site in l
            lop[i, i] += lf(l, site)
            i += 1
        end
    catch e
        if e isa DimensionMismatch
            error("check lambda return type; should be AbstractMatrix")
        else
            rethrow()
        end
    end
    lop
end
function _diag_operator!(lop::LatticeOperator, m::AbstractMatrix)
    for i in 1:length(lop.basis.lattice)
        lop[i, i] += m
    end
    lop
end

diag_operator(f::Function, l::Lattice) =
    _diag_operator!(_zero_on_basis(l, f), _propagate_lattice_args(f, l))
diag_operator(l::Lattice, m::AbstractMatrix) =
    _diag_operator!(_zero_on_basis(l, m), m)
function diag_operator(f::Function, bas::Basis)
    N = bas.internal_dim
    one = Matrix(I, N, N)
    lf = _propagate_lattice_args(f, bas.lattice)
    _diag_operator!(
        _zero_on_basis(bas), (l, site) -> lf(l, site)::Number * one)
end

function coord_operators(b::Basis)
    N = b.internal_dim
    d = dims(b.lattice)
    i = 1
    eye = Matrix(I, (N, N))
    xyz_operators = [LatticeVecOrMat(b, op_mat)
        for op_mat
        in eachslice(zeros(length(b), length(b), d), dims=3)]
    for site in b.lattice
        crd = coords(b.lattice, site)
        for j in 1:d
            xyz_operators[j][i, i] = crd[j] * eye
        end
        i += 1
    end
    return xyz_operators
end

coord_operators(l::Lattice, N::Int) = coord_operators(Basis(l, N))

diag_aggregate(f::Function, lo::LatticeVecOrMat) =
    LatticeValue(lo.basis.lattice, [f(lo[i, i]) for i in 1:length(lo.basis.lattice)])

@inline _make_wrapper(op, ::Nothing) = op
@inline _make_wrapper(op::Number, ::Basis) = op
@inline _make_wrapper(op::Any, b::Basis) = LatticeVecOrMat(b, op)

@inline _unwrap_from_macro(f, args...; kw...) = _unwrap(f, args; kw...)
@inline _unwrap(T::Type, args::Tuple; kw...) = T(args...; kw...)

@inline _unwrap(f::Function, args::Tuple; kw...) = _unwrap(f, (), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap(f, (checked_args..., el), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap(f::Function, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...) =
    _unwrap_wlattice(f, el.basis, (checked_args..., el.operator), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)

@inline _unwrap_nolattice(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)
@inline function _unwrap_nolattice(f::Function, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...)
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
@inline function _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...)
    el.basis != basis && throw(ArgumentError("basis mismatch:\n$(el.basis)\n$basis"))
    _unwrap_wlattice(f, basis, (checked_args..., el.operator), args; kw...)
end
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, checked_args, args[begin], Base.tail(args); kw...)
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, ::Tuple{}; kw...) =
    _make_wrapper(f(checked_args...; kw...), basis)

LatticeSummable = Union{LatticeVecOrMat,UniformScaling}
import Base: +, -, *, /, ^, adjoint, copy
@inline +(los::LatticeSummable...) = _unwrap(+, los)
@inline -(lo::LatticeVecOrMat) = _unwrap(-, (lo,))
@inline -(lo1::LatticeSummable, lo2::LatticeSummable) = _unwrap(-, (lo1, lo2))
@inline *(los::LatticeVecOrMat...) = _unwrap(*, los)
@inline *(num::Number, lo::LatticeVecOrMat) = _unwrap(*, (num, lo))
@inline *(lo::LatticeVecOrMat, num::Number) = _unwrap(*, (lo, num))
@inline /(lo::LatticeVecOrMat, num::Number) = _unwrap(/, (lo, num))
@inline ^(lo::LatticeVecOrMat, num::Number) = _unwrap(^, (lo, num))
@inline adjoint(lo::LatticeVecOrMat) = _unwrap(adjoint, (lo,))
@inline copy(lo::LatticeVecOrMat) = _unwrap(copy, (lo,))

@inline _get_internal(ao::AbstractVecOrMat) = ao
@inline _get_basis(_::AbstractVecOrMat) = nothing
@inline _get_internal(lo::LatticeVecOrMat) = lo.operator
@inline _get_basis(lo::LatticeVecOrMat) = lo.basis

_wrap_smart!(expr::Any) = expr
function _wrap_smart!(expr::Expr)
    Meta.isexpr(expr, :escape) && error("do not use @on_lattice macro in other macros")
    local _begin = 1
    if Meta.isexpr(expr, :call)
        if expr.args[1] === :.|>
            lat_sym, fn_sym = expr.args[2:3]
            expr.args = Any[
                :_make_wrapper,
                :(_get_internal($(esc(lat_sym))) .|> $(esc(fn_sym))),
                :(_get_basis($(esc(lat_sym))))
            ]
            return expr
        end
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
