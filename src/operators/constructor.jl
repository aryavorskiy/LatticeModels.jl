import QuantumOpticsBase: Operator, check_samebases

struct MatrixWrite{T}
    val::T
    overwrite::Bool
end
Base.:(*)(mw::MatrixWrite, n::Number) = MatrixWrite(mw.val * n, mw.overw)
function combine_writes(mw1::MatrixWrite, mw2::MatrixWrite)
    if mw2.overwrite
        return mw2
    else
        return MatrixWrite(mw1.val + mw2.val, true)
    end
end
Base.convert(::Type{MatrixWrite{T}}, mw::MatrixWrite) where T = MatrixWrite(convert(T, mw.val), mw.overwrite)
Base.zero(::Type{MatrixWrite{T}}) where T = MatrixWrite(zero(T), true)
to_number(mw::MatrixWrite) = mw.val

struct SparseMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{MatrixWrite{T}}
    SparseMatrixBuilder{T}(x::Int, y::Int) where T = new{T}((x, y), Int[], Int[], MatrixWrite{T}[])
end
function to_matrix(bh::SparseMatrixBuilder)
    _mat = sparse(bh.Is, bh.Js, bh.Vs, bh.size..., combine_writes)
    return SparseMatrixCSC(_mat.m, _mat.n, _mat.colptr, _mat.rowval, to_number.(_mat.nzval))
end

Base.@propagate_inbounds function Base.setindex!(bh::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int; overwrite=true, factor=1)
    # number increment
    push!(bh.Is, i1)
    push!(bh.Js, i2)
    push!(bh.Vs, MatrixWrite(rhs * factor, overwrite))
    return nothing
end
Base.@propagate_inbounds function Base.setindex!(bh::SparseMatrixBuilder, rhs::AbstractMatrix,
        is1::AbstractArray, is2::AbstractArray; kw...)
    for (i1, j1) in enumerate(is1), (i2, j2) in enumerate(is2)
        v = rhs[i1, i2]
        iszero(v) || Base.setindex!(bh, v, j1, j2; kw...)
    end
end
Base.@propagate_inbounds function increment!(bh::SparseMatrixBuilder, rhs::SparseMatrixCSC)
    # global increment
    nis, njs, nvs = findnz(rhs)
    append!(bh.Is, nis)
    append!(bh.Js, njs)
    append!(bh.Vs, MatrixWrite.(nvs * factor, false))
end

struct UniversalView{AT}
    axes::AT
end
struct WrappedUniversalView{RT,AT,NT}
    axes::AT
    rhs::RT
    factor::NT
end
Base.:(+)(uv::UniversalView, rhs) = WrappedUniversalView(uv.axes, rhs, 1)
Base.:(-)(uv::UniversalView, rhs) = WrappedUniversalView(uv.axes, rhs, -1)
Base.adjoint(uvw::WrappedUniversalView) = WrappedUniversalView(reverse(uvw.axes), uvw.rhs', uvw.factor')

Base.getindex(::SparseMatrixBuilder, I...) = UniversalView(I)
Base.@propagate_inbounds function Base.setindex!(bh::SparseMatrixBuilder, wuv::WrappedUniversalView, I...; factor=1)
    @boundscheck wuv.axes != I && throw(ArgumentError("Invalid indexing. Only `A[I...] (+=|-=|=) B` is allowed."))
    Base.setindex!(bh, wuv.rhs, I...; overwrite=false, factor=wuv.factor * factor)
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
    @boundscheck all(==(internal_length(sample)), size(mat))
    return mat
end
op_to_matrix(::Sample, op) =
    throw(ArgumentError("unsupported on-site operator type $(typeof(op))"))

# convert site1=>site2 pair to (i,j,fact) tuple
function expand_bond(l::AbstractLattice, site1::AbstractSite, site2::AbstractSite, field::AbstractField)
    ifact, new_site1 = shift_site(l, site1)
    jfact, new_site2 = shift_site(l, site2)
    i = site_index(l, new_site1)
    j = site_index(l, new_site2)
    i === nothing && return nothing
    j === nothing && return nothing
    field_fact = exp(-2π * im * line_integral(field, site1.coords, site2.coords))
    !isfinite(field_fact) && @warn("got NaN or Inf when finding the phase factor")
    return i, j, jfact * field_fact * ifact'
end

Base.getindex(::OperatorBuilder, axes...) = UniversalView(axes)

function _preprocess_op(sample::Sample, uvw::WrappedUniversalView, axes, new_axes)
    new_mat = op_to_matrix(sample, uvw.rhs)
    @boundscheck if uvw.axes != axes
        throw(ArgumentError("Invalid `OperatorBuilder` indexing. Only `A[site1, site2] (+=|-=|=) B` allowed, where B is a number, matrix or `Operator`"))
    end
    return WrappedUniversalView(new_axes, new_mat, uvw.factor)
end
_preprocess_op(sample::Sample, op, _, _) = op_to_matrix(sample, op)
Base.@propagate_inbounds function Base.setindex!(opbuilder::OperatorBuilder, rhs, site1::AbstractSite, site2::AbstractSite)
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
