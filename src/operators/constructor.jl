import QuantumOpticsBase: Operator, check_samebases

struct ArrayEntry{T}
    val::T
    overwrite::Bool
end
Base.:(*)(mw::ArrayEntry, n::Number) = ArrayEntry(mw.val * n, mw.overw)
function combine_writes(mw1::ArrayEntry, mw2::ArrayEntry)
    if mw2.overwrite
        return mw2
    else
        return ArrayEntry(mw1.val + mw2.val, true)
    end
end
Base.convert(::Type{ArrayEntry{T}}, mw::ArrayEntry) where T = ArrayEntry(convert(T, mw.val), mw.overwrite)
Base.zero(::Type{ArrayEntry{T}}) where T = ArrayEntry(zero(T), true)
to_number(mw::ArrayEntry) = mw.val

struct SparseMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{ArrayEntry{T}}
    SparseMatrixBuilder{T}(x::Int, y::Int) where T = new{T}((x, y), Int[], Int[], ArrayEntry{T}[])
end
function to_matrix(A::SparseMatrixBuilder)
    _mat = sparse(A.Is, A.Js, A.Vs, A.size..., combine_writes)
    return SparseMatrixCSC(_mat.m, _mat.n, _mat.colptr, _mat.rowval, to_number.(_mat.nzval))
end

Base.@propagate_inbounds function Base.setindex!(A::SparseMatrixBuilder, B::Number, i1::Int, i2::Int; overwrite=true, factor=1)
    # number increment
    push!(A.Is, i1)
    push!(A.Js, i2)
    push!(A.Vs, ArrayEntry(B * factor, overwrite))
    return nothing
end
Base.@propagate_inbounds function Base.setindex!(A::SparseMatrixBuilder, B::AbstractMatrix,
        is1::AbstractArray, is2::AbstractArray; overwrite=true, factor=1)
    for (i1, j1) in enumerate(is1), (i2, j2) in enumerate(is2)
        v = B[i1, i2]
        iszero(v) || (A[j1, j2, overwrite=overwrite, factor=factor] = v)
    end
end
Base.@propagate_inbounds function increment!(A::SparseMatrixBuilder, B::SparseMatrixCSC)
    # global increment
    nis, njs, nvs = findnz(B)
    append!(A.Is, nis)
    append!(A.Js, njs)
    append!(A.Vs, ArrayEntry.(nvs * factor, false))
end

struct BuilderView{AT}
    axes::AT
end
struct BuilderIncrementIndex{RT,AT,NT}
    axes::AT
    B::RT
    factor::NT
end
Base.:(+)(v::BuilderView, B) = BuilderIncrementIndex(v.axes, B, 1)
Base.:(-)(v::BuilderView, B) = BuilderIncrementIndex(v.axes, B, -1)
Base.adjoint(incr::BuilderIncrementIndex) = BuilderIncrementIndex(reverse(incr.axes), incr.B', incr.factor')

Base.getindex(::SparseMatrixBuilder, I...) = BuilderView(I)
Base.@propagate_inbounds function Base.setindex!(b::SparseMatrixBuilder, incr::BuilderIncrementIndex, I...; factor=1)
    @boundscheck incr.axes != I && throw(ArgumentError("Invalid indexing. Only `A[I...] (+=|-=|=) B` is allowed."))
    b[I..., overwrite=false, factor=incr.factor * factor] = incr.B
end

struct OperatorBuilder{SystemT, FieldT, T}
    sys::SystemT
    field::FieldT
    internal_length::Int
    mat_builder::SparseMatrixBuilder{T}
    auto_hermitian::Bool
end

function OperatorBuilder(T::Type{<:Number}, sys::SystemT;
        field::FieldT=NoField(), auto_hermitian=false) where {SystemT<:System, FieldT<:AbstractField}
    oneparticle_len = length(onebodybasis(sys))
    OperatorBuilder{SystemT, FieldT, T}(sys, field, internal_length(sys),
        SparseMatrixBuilder{T}(oneparticle_len, oneparticle_len), auto_hermitian)
end
@accepts_system_t OperatorBuilder

sample(opb::OperatorBuilder) = sample(opb.sys)
const OpBuilderWithInternal = OperatorBuilder{<:System{<:SampleWithInternal}}
const OpBuilderWithoutInternal = OperatorBuilder{<:System{<:SampleWithoutInternal}}

_internal_one_mat(sample::SampleWithInternal) = internal_one(sample).data
_internal_one_mat(::SampleWithoutInternal) = SMatrix{1,1}(1)
op_to_matrix(sample::SampleWithoutInternal, n::Number) = n * _internal_one_mat(sample)
function op_to_matrix(sample::SampleWithInternal, op::DataOperator)
    @boundscheck QuantumOpticsBase.check_samebases(basis(op), internal_basis(sample))
    return op.data
end
function op_to_matrix(sample::Sample, mat::AbstractMatrix)
    @boundscheck all(==(internal_length(sample)), size(mat)) ||
        throw(ArgumentError("matrix size does not match on-site dims"))
    return mat
end
op_to_matrix(::Sample, op) =
    throw(ArgumentError("cannot interpret $(typeof(op)) as on-site operator"))

function expand_bond(l::AbstractLattice, site1::AbstractSite, site2::AbstractSite, field::AbstractField)
    s1 = resolve_site(l, site1)
    s2 = resolve_site(l, site2)
    s1 === nothing && return nothing
    s2 === nothing && return nothing
    return expand_bond(l, s1, s2, field)
end

@inline function expand_bond(::AbstractLattice, s1::ResolvedSite, s2::ResolvedSite, field::AbstractField)
    field_fact = exp(-2π * im * line_integral(field, s1.site.coords, s2.site.coords))
    return s1.index, s2.index, s2.factor * field_fact * s1.factor'
end

Base.getindex(::OperatorBuilder, axes...) = BuilderView(axes)

function _preprocess_op(sample::Sample, uvw::BuilderIncrementIndex, axes, new_axes)
    new_mat = op_to_matrix(sample, uvw.B)
    @boundscheck if uvw.axes != axes
        throw(ArgumentError("Invalid `OperatorBuilder` indexing. Only `A[site1, site2] (+=|-=|=) B` allowed, where B is a number, matrix or `Operator`"))
    end
    return BuilderIncrementIndex(new_axes, new_mat, uvw.factor)
end
_preprocess_op(sample::Sample, op, _, _) = op_to_matrix(sample, op)
Base.@propagate_inbounds function Base.setindex!(opbuilder::OperatorBuilder, rhs, site1, site2)
    ijfact = expand_bond(lattice(opbuilder), site1, site2, opbuilder.field)
    ijfact === nothing && return nothing
    i, j, total_factor = ijfact
    N = internal_length(opbuilder)
    is = (i - 1) * N+1:i * N
    js = (j - 1) * N+1:j * N
    new_rhs = _preprocess_op(sample(opbuilder), rhs, (site1, site2), (is, js))
    opbuilder.mat_builder[is, js, factor = total_factor] = new_rhs
    if opbuilder.auto_hermitian && i != j
        opbuilder.mat_builder[js, is, factor = total_factor'] = new_rhs'
    end
    return nothing
end

_build_manybody_maybe(sys::ManyBodySystem, op::AbstractOperator) = manybodyoperator(sys, op)
_build_manybody_maybe(::OneParticleBasisSystem, op::AbstractOperator) = op
function QuantumOpticsBase.Operator(opb::OperatorBuilder; warning=true)
    op = Operator(onebodybasis(opb.sys), to_matrix(opb.mat_builder))
    if warning && !opb.auto_hermitian && !ishermitian(op)
        @warn "The resulting operator is not hermitian. Set `warning=false` to hide this message or add `auto_hermitian=true` to the `OperatorBuilder` constructor."
    end
    return _build_manybody_maybe(opb.sys, op)
end
Hamiltonian(opb::OperatorBuilder; kw...) = Hamiltonian(opb.sys, Operator(opb; kw...))
