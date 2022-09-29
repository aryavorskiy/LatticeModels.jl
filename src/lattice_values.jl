using LinearAlgebra, Statistics, Logging
import Base: length, iterate, show

struct Basis
    lattice::AbstractLattice
    internal_dim::Int
end

length(b::Basis) = length(b.lattice) * b.internal_dim
dims_internal(b::Basis) = b.internal_dim

function show(io::IO, m::MIME"text/plain", b::Basis)
    println(io, "Basis with $(b.internal_dim)-dimensional internal phase space")
    print(io, "on ")
    show(io, m, b.lattice)
end

struct LatticeValue{VT<:AbstractVector}
    lattice::AbstractLattice
    vector::VT
end

LatticeValue(f, l::AbstractLattice) =
    (lf = _propagated_lattice_args(f, l); LatticeValue(l, [lf(l, site) for site in l]))

LatticeValue(f, VT::Type, l::AbstractLattice) =
    (lf = _propagated_lattice_args(f, l);
    LatticeValue(l, convert(VT, [lf(l, site) for site in l])))

function show(io::IO, m::MIME"text/plain", lv::LatticeValue{VT}) where VT
    println(io, "LatticeValue with inner type $VT")
    print(io, "on ")
    show(io, m, lv.lattice)
end

convert_inner_type(VT::Type, lv::LatticeValue) =
    LatticeValue(lv.lattice, convert(VT, lv.vector))

length(lv::LatticeValue) = length(lv.vector)
iterate(lv::LatticeValue, args...) = iterate(lv.vector, args...)

@recipe function f(lv::LatticeValue)
    lv.lattice, lv.vector
end

struct LatticeVecOrMat{MT<:AbstractVecOrMat}
    operator::MT
    basis::Basis
    function LatticeVecOrMat(basis::Basis, operator::AbstractVecOrMat)
        @assert all(size(operator) .== length(basis))
        new{typeof(operator)}(operator, basis)
    end
end

LatticeOperator{T} = LatticeVecOrMat{T} where T<: AbstractMatrix

size(lv::LatticeVecOrMat) = size(lv.operator)
dims_internal(lv::LatticeVecOrMat) = lv.basis.internal_dim

function show(io::IO, m::MIME"text/plain", lv::LatticeVecOrMat{MT}) where MT<:AbstractMatrix
    println(io, join(size(lv), "×") * " LatticeMatrix with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeVecOrMat{MT}) where MT<:AbstractVector
    println(io, "$(length(lv.operator))-length LatticeVector with inner type $MT")
    print(io, "on ")
    show(io, m, lv.basis)
end

convert_inner_type(MT::Type, lo::LatticeVecOrMat) =
    LatticeValue(lo.basis, convert(MT, lo.operator))

function diag_operator(f::Function, l::AbstractLattice)
    ops = _map_to_lattice(f, l)
    @assert allequal(size(op) for op in ops) "inconsistent size of lambda return value"
    if size(ops[1]) == ()
        throw(ArgumentError("lambda returns a number, not a matrix"))
    end
    N = size(ops[1])[begin]
    # matrix = zeros(ComplexF64, (N * length(l), N * length(l)))
    matrix = zero(similar(ops[1], ComplexF64, (N * length(l), N * length(l))))
    for i in 1:length(l)
        matrix[N * (i - 1) + 1: N * i, N * (i - 1) + 1: N * i] .= ops[i]
    end
    LatticeVecOrMat(Basis(l, N), matrix)
end

function diag_operator(l::AbstractLattice, m::Matrix)
    N = size(m)[begin]
    matrix = zero(similar(m, ComplexF64, (N * length(l), N * length(l))))
    for i in 1:length(l)
        matrix[N * (i - 1) + 1: N * i, N * (i - 1) + 1: N * i] .= m
    end
    LatticeVecOrMat(Basis(l, N), matrix)
end
diag_operator(n::Number, l::AbstractLattice) = diag_operator([n;;], l)

function coord_operators(b::Basis)
    N = b.internal_dim
    d = dims(b.lattice)
    i = 1
    eye = Matrix(I, (N, N))
    xyz_operators = zeros(length(b), length(b), d)
    for site in b.lattice
        crd = coords(b.lattice, site)
        for j in 1:d
            xyz_operators[N * (i - 1) + 1: N * i, N * (i - 1) + 1: N * i, j] = crd[j] * eye
        end
        i += 1
    end
    return [LatticeVecOrMat(b, op_mat) for op_mat in eachslice(xyz_operators, dims=3)]
end

coord_operators(l::AbstractLattice, N::Int) = coord_operators(Basis(l, N))

function diag_aggregate(f::Function, lo::LatticeVecOrMat)
    local N = lo.basis.internal_dim
    LatticeValue(
        lo.basis.lattice,
        [f(lo.operator[N * (i - 1) + 1: N * i, N * (i - 1) + 1: N * i])
            for i in 1:length(lo.basis.lattice)]
    )
end

@inline _get_internal(ao::AbstractVecOrMat) = ao
@inline _get_basis(_::AbstractVecOrMat) = nothing
@inline _get_internal(lo::LatticeVecOrMat) = lo.operator
@inline _get_basis(lo::LatticeVecOrMat) = lo.basis
@inline _get_internal(lv::LatticeValue) = lv.vector
@inline _get_basis(lv::LatticeValue) = lv.lattice

@inline _make_wrapper(op, ::Nothing) = op
@inline _make_wrapper(op::Number, ::Basis) = op
@inline _make_wrapper(op::Any, b::Basis) = LatticeVecOrMat(b, op)
@inline _make_wrapper(op, l::AbstractLattice) = LatticeValue(l, op)

@inline _unwrap_from_macro(f, args...; kw...) = _unwrap(f, args; kw...)
@inline _unwrap(T::Type, args::Tuple; kw...) = T(args...; kw...)

@inline _unwrap(f::Function, args::Tuple; kw...) = _unwrap(f, (), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap(f, (checked_args..., el), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap(f, checked_args, args[1], Base.tail(args); kw...)
@inline _unwrap(f::Function, checked_args::Tuple, ::Tuple{}; kw...) = f(checked_args...; kw...)

@inline _unwrap(f::Function, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...) =
    _unwrap_wlattice(f, _get_basis(el), (checked_args..., _get_internal(el)), args; kw...)
@inline _unwrap(f::Function, checked_args::Tuple, el::AbstractVecOrMat, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)

@inline _unwrap_nolattice(f::Function, checked_args::Tuple, el::Any, args::Tuple; kw...) =
    _unwrap_nolattice(f, (checked_args..., el), args; kw...)
@inline function _unwrap_nolattice(f::Function, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...)
    @warn "avoid using lattice operators and matrices in one function call"
    _unwrap_wlattice(f, _get_basis(el), (checked_args..., _get_internal(el)), args; kw...)
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
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, el::LatticeVecOrMat, args::Tuple; kw...) =
    _get_basis(el) != basis ?
    throw(ArgumentError(
        "basis does not match one of operands before:
        $(_get_basis(el))
        $basis")) :
    _unwrap_wlattice(f, basis, (checked_args..., _get_internal(el)), args; kw...)
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, args::Tuple; kw...) =
    _unwrap_wlattice(f, basis, checked_args, args[begin], Base.tail(args); kw...)
@inline _unwrap_wlattice(f::Function, basis::Any, checked_args::Tuple, ::Tuple{}; kw...) =
    _make_wrapper(f(checked_args...; kw...), basis)

LatticeSummable = Union{LatticeVecOrMat, UniformScaling}
import Base: +, -, *, /, adjoint, copy
@inline +(los::LatticeSummable...) = _unwrap(+, los)
@inline -(lo::LatticeVecOrMat) = _unwrap(-, (lo,))
@inline -(lo1::LatticeSummable, lo2::LatticeSummable) = _unwrap(-, (lo1, lo2))
@inline *(los::LatticeVecOrMat...) = _unwrap(*, los)
@inline *(num::Number, lo::LatticeVecOrMat) = _unwrap(*, (num, lo))
@inline *(lo::LatticeVecOrMat, num::Number) = _unwrap(*, (lo, num))
@inline /(lo::LatticeVecOrMat, num::Number) = _unwrap(/, (lo, num))
@inline adjoint(lo::LatticeVecOrMat) = _unwrap(adjoint, (lo,))
@inline copy(lo::LatticeVecOrMat) = _unwrap(copy, (lo,))

_wrap_smart!(expr::Any) = expr
function _wrap_smart!(expr::Expr)
    Meta.isexpr(expr, :escape) && error("do not use @on_lattice macro in other macros")
    local _begin = 1
    if Meta.isexpr(expr, :call)
        if expr.args[1] === :.|>
            lat_sym, fn_sym =  expr.args[2:3]
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
