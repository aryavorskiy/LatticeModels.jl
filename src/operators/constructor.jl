import QuantumOpticsBase: Operator, check_samebases

Base.@propagate_inbounds function increment!(lhs, rhs, arg1, args...)
    lhs[arg1, args...] += rhs
    return lhs
end
Base.@propagate_inbounds increment!(lhs, rhs) = lhs += rhs; return lhs
Base.@propagate_inbounds function increment!(builder::AbstractArray, rhs::Number, i1::Int, i2::Int; factor=1)
    builder[i1, i2] += rhs * factor
    return nothing
end
Base.@propagate_inbounds function increment!(builder::AbstractArray, rhs::SparseMatrixCSC; factor=1)
    @. builder += rhs * factor
    return builder
end

function _process_increment(expr::Expr)
    !Meta.isexpr(expr, :(+=)) && error("+= increment expected")
    lhs, rhs = expr.args
    if Meta.isexpr(lhs, :ref)
        builder, refs = lhs.args[1], lhs.args[2:end]
        return :(LatticeModels.increment!($builder, $rhs, $(refs...)))
    elseif lhs isa Symbol
        return :($lhs = LatticeModels.increment!($lhs, $rhs))
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
    SparseMatrixBuilder{T}(x::Int, y::Int=x) where T = new{T}((x, y), Int[], Int[], T[])
end
Base.size(smb::SparseMatrixBuilder) = smb.size
to_matrix(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)
to_matrix(mat::AbstractMatrix) = mat

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int; factor=1)
    # number increment
    push!(builder.Is, i1)
    push!(builder.Js, i2)
    push!(builder.Vs, rhs * factor)
    return nothing
end
Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::SparseMatrixCSC; factor=1)
    # global increment
    nis, njs, nvs = findnz(rhs)
    append!(builder.Is, nis)
    append!(builder.Js, njs)
    append!(builder.Vs, nvs * factor)
    return builder
end
Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::AbstractMatrix, is1, is2; kw...)
    for (i1, j1) in enumerate(is1), (i2, j2) in enumerate(is2)
        v = rhs[i1, i2]
        iszero(v) || increment!(builder, v, j1, j2; kw...)
    end
    return builder
end

struct OperatorBuilder{SystemT, FieldT, BT}
    sys::SystemT
    field::FieldT
    internal_length::Int
    builder::BT
    auto_hermitian::Bool
end

function OperatorBuilder(T::Type{<:Number}, sys::SystemT;
        field::FieldT=NoField(), auto_hermitian=false) where {SystemT<:System, FieldT<:AbstractField}
    mat = zeros(T, length(onebodybasis(sys)), length(onebodybasis(sys)))
    OperatorBuilder{SystemT, FieldT, typeof(mat)}(sys, field, internal_length(sys), mat, auto_hermitian)
end
@accepts_system_t OperatorBuilder

const SparseOperatorBuilder{SystemT, FieldT, T} =
    OperatorBuilder{SystemT, FieldT, <:SparseMatrixCSC{T}}
function SparseOperatorBuilder(T::Type{<:Number}, sys::SystemT; field::FieldT=NoField(), auto_hermitian=false
    ) where {SystemT<:System, FieldT<:AbstractField}
    mat = spzeros(T, length(onebodybasis(sys)), length(onebodybasis(sys)))
    OperatorBuilder{SystemT, FieldT, typeof(mat)}(sys, field, internal_length(sys), mat, auto_hermitian)
end
@accepts_system_t SparseOperatorBuilder

const FastSparseOperatorBuilder{SystemT, FieldT, T} =
    OperatorBuilder{SystemT, FieldT, SparseMatrixBuilder{T}}
function FastSparseOperatorBuilder(T::Type{<:Number}, sys::SystemT; field::FieldT=NoField(),
    auto_hermitian=false) where {SystemT<:System, FieldT<:AbstractField}
    mat = SparseMatrixBuilder{T}(length(onebodybasis(sys)))
    OperatorBuilder{SystemT, FieldT, typeof(mat)}(sys, field, internal_length(sys), mat, auto_hermitian)
end
@accepts_system_t FastSparseOperatorBuilder
sample(opb::OperatorBuilder) = sample(opb.sys)
const OpBuilderWithInternal = OperatorBuilder{<:System{<:SampleWithInternal}}
const OpBuilderWithoutInternal = OperatorBuilder{<:System{<:SampleWithoutInternal}}

struct NoMatrixElement end
Base.:(+)(::NoMatrixElement, _) = NoMatrixElement()
Base.:(*)(::NoMatrixElement, _) = NoMatrixElement()

_internal_one_mat(sample::SampleWithInternal) = internal_one(sample).data
_internal_one_mat(sample::SampleWithoutInternal) = 1
op_to_matrix(sample::SampleWithoutInternal, n::Number) = n * _internal_one_mat(sample)
function op_to_matrix(sample::SampleWithInternal, op::DataOperator)
    @boundscheck QuantumOpticsBase.check_samebases(basis(op), internal_basis(sample))
    return op.data
end
op_to_matrix(::Sample, mat::AbstractMatrix) = mat
op_to_matrix(::Sample, op) =
    throw(ArgumentError("unsupported on-site operator type $(typeof(op))"))

# convert site1=>site2 pair to (i,j,fact) tuple
function expand_bond(l::AbstractLattice, bond::SingleBond, field::AbstractField)
    site1, site2 = bond
    ifact, new_site1 = shift_site(l, site1)
    jfact, new_site2 = shift_site(l, site2)
    i = site_index(l, new_site1)
    j = site_index(l, new_site2)
    i === nothing && return nothing
    j === nothing && return nothing
    field_fact = exp(-2π * im * line_integral(field, site1.coords, site2.coords))
    return i, j, jfact * field_fact * ifact'
end

Base.@propagate_inbounds function increment!(opbuilder::OperatorBuilder, rhs, site1::AbstractSite, site2::AbstractSite)
    new_rhs = op_to_matrix(sample(opbuilder), rhs)
    ijfact = expand_bond(lattice(opbuilder), site1 => site2, opbuilder.field)
    ijfact === nothing && return nothing
    N = internal_length(opbuilder)
    increment!(opbuilder.builder, new_rhs, (i-1)*N+1:i*N, (j-1)*N+1:j*N, factor = total_factor)
    if opbuilder.auto_hermitian && i != j
        increment!(opbuilder.builder, new_rhs', (j-1)*N+1:j*N, (i-1)*N+1:i*N, factor = total_factor')
    end
    return nothing
end

Base.@propagate_inbounds function Base.getindex(opbuilder::OperatorBuilder, site1::AbstractSite, site2::AbstractSite)
    ijfact = expand_bond(lattice(opbuilder), site1 => site2, opbuilder.field)
    ijfact === nothing && return NoMatrixElement()
    i, j, total_factor = ijfact
    N = internal_length(opbuilder)
    mat = opbuilder.builder[(i-1)*N+1:i*N, (j-1)*N+1:j*N] / total_factor
    if hasinternal(opbuilder)
        return Operator(internal_basis(opbuilder), mat)
    else
        return only(mat)
    end
end
Base.getindex(::FastSparseOperatorBuilder, ::AbstractSite, ::AbstractSite) =
    error("`FastSparseOperatorBuilder` does not support normal indexing. Maybe you forgot to use `@increment` macro?")

Base.@propagate_inbounds function Base.setindex!(opbuilder::OperatorBuilder, rhs, site1::AbstractSite, site2::AbstractSite)
    new_rhs = op_to_matrix(sample(opbuilder), rhs)
    ijfact = expand_bond(lattice(opbuilder), site1 => site2, opbuilder.field)
    ijfact === nothing && return NoMatrixElement()
    i, j, total_factor = ijfact
    N = internal_length(opbuilder)
    @boundscheck !all(size(new_rhs) == (N, N)) &&
        error("rhs size mismatch: Expected $N×$N, got $(join("×", size(new_rhs))))")
    opbuilder.builder[(i-1)*N+1:i*N, (j-1)*N+1:j*N] = new_rhs * total_factor
    if opbuilder.auto_hermitian && i != j
        opbuilder.builder[(j-1)*N+1:j*N, (i-1)*N+1:i*N] = new_rhs' * total_factor'
    end
    return rhs
end
Base.setindex!(::FastSparseOperatorBuilder, ::Any, ::AbstractSite, ::AbstractSite) =
    error("`FastSparseOperatorBuilder` does not support normal indexing. Use increments `+=` and `@increment` macro.")
Base.setindex!(::OperatorBuilder, ::NoMatrixElement, ::AbstractSite, ::AbstractSite) = NoMatrixElement()

function _build_manybody_maybe(ps::NParticles, op::AbstractOperator)
    check_samebases(onebodybasis(ps), basis(op))
    return manybodyoperator(ManyBodyBasis(basis(op), occupations(ps)), op)
end
_build_manybody_maybe(::OneParticleBasisSystem, op::AbstractOperator) = op
function QuantumOpticsBase.Operator(opb::OperatorBuilder; warning=true)
    op = Operator(onebodybasis(opb.sys), to_matrix(opb.builder))
    if warning && !opb.auto_hermitian && !ishermitian(op)
        @warn "The resulting operator is not hermitian. Set `warning=false` to hide this message or add `auto_hermitian=true` to the `OperatorBuilder` constructor."
    end
    return _build_manybody_maybe(opb.sys, op)
end
Hamiltonian(opb::OperatorBuilder; kw...) = Hamiltonian(opb.sys, Operator(opb; kw...))
