import QuantumOpticsBase: Operator, check_samebases

Base.@propagate_inbounds function increment!(builder::AbstractArray, rhs::Number, i1::Int, i2::Int; factor=1)
    builder[i1, i2] += rhs * factor
    return nothing
end
Base.@propagate_inbounds function increment!(builder::AbstractArray, rhs::SparseMatrixCSC; factor=1)
    @. builder += rhs * factor
    return nothing
end

function _process_increment(expr::Expr)
    !Meta.isexpr(expr, :(+=)) && error("+= increment expected")
    lhs, rhs = expr.args
    if Meta.isexpr(lhs, :ref)
        builder, refs = lhs.args[1], lhs.args[2:end]
        return :(LatticeModels.increment!($builder, $rhs, $(refs...)))
    elseif lhs isa Symbol
        return :(LatticeModels.increment!($lhs, $rhs))
    end
end
function _process_increment_recursive!(expr)
    !Meta.isexpr(expr, (:if, :elseif, :for, :while, :block)) && return expr
    for i in eachindex(expr.args)
        if Meta.isexpr(expr.args[i], :(+=))
            expr.args[i] = _process_increment(expr.args[i])
        else
            _process_increment_recursive!(expr.args[i])
        end
    end
    return expr
end
macro increment(expr)
    if Meta.isexpr(expr, :(+=))
        return esc(_process_increment(expr))
    else
        return esc(_process_increment_recursive!(expr))
    end
end

struct SparseMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    SparseMatrixBuilder{T}(sz) where T = new{T}(sz, Int[], Int[], T[])
end
SparseMatrixBuilder{T}(sz...) where T = SparseMatrixBuilder{T}(sz)
SparseMatrixBuilder{T}(size_x::Int) where T = SparseMatrixBuilder{T}(size_x, size_x)
Base.size(smb::SparseMatrixBuilder) = smb.size
to_matrix(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)
to_matrix(mat::AbstractMatrix) = mat

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int; factor=1)
    push!(builder.Is, i1)
    push!(builder.Js, i2)
    push!(builder.Vs, rhs * factor)
    return rhs
end
Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::SparseMatrixCSC; factor=1)
    nis, njs, nvs = findnz(rhs)
    append!(builder.Is, nis)
    append!(builder.Js, njs)
    append!(builder.Vs, nvs * factor)
    return rhs
end

const AbstractMatrixBuilder = Union{AbstractMatrix, SparseMatrixBuilder}
Base.@propagate_inbounds function increment!(builder::AbstractMatrixBuilder, rhs::AbstractMatrix, i1::Int, i2::Int; kw...)
    N = size(rhs)[1]
    for i in 1:N, j in 1:N
        v = rhs[i, j]
        iszero(v) || increment!(builder, v, i + N * (i1 - 1), j + N * (i2 - 1); kw...)
    end
    return rhs
end
Base.@propagate_inbounds increment!(builder::AbstractMatrixBuilder, rhs::DataOperator, i1::Int, i2::Int; kw...) =
    increment!(builder, rhs.data, i1, i2; kw...)

struct OperatorBuilder{SystemT, FieldT, BT}
    sys::SystemT
    field::FieldT
    internal_length::Int
    builder::BT
    auto_hermitian::Bool
end

function OperatorBuilder(T::Type{<:Number}, sys::SystemT; field::FieldT=NoField(), auto_hermitian=false
        ) where {SystemT<:System, FieldT<:AbstractField}
    mat = zeros(T, length(onebodybasis(sys)), length(onebodybasis(sys)))
    OperatorBuilder{SystemT, FieldT, typeof(mat)}(sys, field, internal_length(sys), mat, auto_hermitian)
end
(::Type{T})(args...; kw...) where {T<:OperatorBuilder} = T(ComplexF64, args...; kw...)
(::Type{T})(::Type, ::Type, args...; kw...) where {T<:OperatorBuilder} =
    throw(MethodError(T, args))
@accepts_system OperatorBuilder

const FastSparseOperatorBuilder{SystemT, FieldT, T} =
    OperatorBuilder{SystemT, FieldT, SparseMatrixBuilder{T}}
function FastSparseOperatorBuilder(T::Type{<:Number}, sys::SystemT; field::FieldT=NoField(),
    auto_hermitian=false) where {SystemT<:System, FieldT<:AbstractField}
    mat = SparseMatrixBuilder{T}(length(onebodybasis(sys)), length(onebodybasis(sys)))
    OperatorBuilder{SystemT, FieldT, typeof(mat)}(sys, field, internal_length(sys), mat, auto_hermitian)
end
@accepts_system FastSparseOperatorBuilder
sample(opb::OperatorBuilder) = sample(opb.sys)
const OpBuilderWithInternal = OperatorBuilder{<:System{<:SampleWithInternal}}
const OpBuilderWithoutInternal = OperatorBuilder{<:System{<:SampleWithoutInternal}}

struct NoMatrixElement end
Base.:(+)(::NoMatrixElement, _) = NoMatrixElement()
Base.:(*)(::NoMatrixElement, _) = NoMatrixElement()

preprocess_rhs(::OpBuilderWithoutInternal, n::Number) = SMatrix{1, 1, Int}(1)
preprocess_rhs(opb::OperatorBuilder, n::Number) = n * internal_one(opb.sys)
function preprocess_rhs(opb::OpBuilderWithInternal, op::DataOperator)
    @boundscheck QuantumOpticsBase.check_samebases(basis(op), internal_basis(opb.sys))
    return op.data
end
preprocess_rhs(::OperatorBuilder, mat::AbstractMatrix) = mat

function preprocess_sites(l::AbstractLattice, field::AbstractField, site1::AbstractSite, site2::AbstractSite)
    ifact, new_site1 = shift_site(l, site1)
    jfact, new_site2 = shift_site(l, site2)
    i = site_index(l, new_site1)
    j = site_index(l, new_site2)
    field_fact = exp(-2π * im * line_integral(field, site1.coords, site2.coords))
    i === nothing && return nothing
    j === nothing && return nothing
    return i, j, jfact * field_fact * ifact'
end

Base.@propagate_inbounds function increment!(opbuilder::OperatorBuilder, rhs, site1, site2)
    new_rhs = preprocess_rhs(opbuilder, rhs)
    prop = preprocess_sites(lattice(opbuilder), opbuilder.field, site1, site2)
    prop === nothing && return NoMatrixElement()
    i, j, total_factor = prop
    increment!(opbuilder.builder, new_rhs, i, j, factor = total_factor)
    if opbuilder.auto_hermitian && i != j
        increment!(opbuilder.builder, new_rhs', j, i, factor = total_factor')
    end
    return rhs
end

Base.@propagate_inbounds function Base.getindex(opbuilder::OperatorBuilder, site1, site2)
    prop = preprocess_sites(lattice(opbuilder), opbuilder.field, site1, site2)
    prop === nothing && return NoMatrixElement()
    i, j, total_factor = prop
    N = internal_length(opbuilder)
    mat = opbuilder.builder[(i - 1) * N + 1:i * N, (j - 1) * N + 1:j * N] * total_factor'
    return Operator(internal_basis(opbuilder), mat)
end

Base.@propagate_inbounds function Base.setindex!(opbuilder::OperatorBuilder, rhs, site1, site2)
    new_rhs = preprocess_rhs(opbuilder, rhs)
    prop = preprocess_sites(lattice(opbuilder), opbuilder.field, site1, site2)
    prop === nothing && return NoMatrixElement()
    i, j, total_factor = prop
    N = internal_length(opbuilder)
    @boundscheck !all(size(new_rhs) == (N, N)) &&
        error("rhs size mismatch: Expected $N×$N, got $(join("×", size(new_rhs))))")
    opbuilder.builder[(i - 1) * N + 1:i * N, (j - 1) * N + 1:j * N] = new_rhs * total_factor
    if opbuilder.auto_hermitian && i != j
        opbuilder.builder[(j - 1) * N + 1:j * N, (i - 1) * N + 1:i * N] = new_rhs' * total_factor'
    end
    return rhs
end
Base.setindex!(::OperatorBuilder, ::NoMatrixElement, site1, site2) = NoMatrixElement()

function QuantumOpticsBase.Operator(opb::OperatorBuilder; warning=true)
    op = Operator(onebodybasis(opb.sys), to_matrix(opb.builder))
    if warning && !opb.auto_hermitian && !ishermitian(op)
        @warn "The resulting operator is not hermitian. Set `warning=false` to hide this message or add `auto_hermitian=true` to the `OperatorBuilder` constructor."
    end
    return manybodyoperator(opb.sys, op)
end
