import QuantumOpticsBase: Operator, manybodyoperator, check_samebases

abstract type AbstractMatrixBuilder{T} end

struct ArrayEntry{T}
    val::T
    overwrite::Bool
end
function combine_entries(mw1::ArrayEntry, mw2::ArrayEntry)
    if mw2.overwrite
        return mw2
    else
        return ArrayEntry(mw1.val + mw2.val, true)
    end
end
Base.convert(::Type{ArrayEntry{T}}, mw::ArrayEntry) where T = ArrayEntry(convert(T, mw.val), mw.overwrite)
Base.convert(::Type{ArrayEntry{T}}, x::Any) where T = ArrayEntry(convert(T, x), false)
Base.zero(::Type{ArrayEntry{T}}) where T = ArrayEntry(zero(T), true)
to_number(mw::ArrayEntry) = mw.val

struct VectorMatrixBuilder{T} <: AbstractMatrixBuilder{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    VectorMatrixBuilder{T}(x::Int, y::Int) where T = new{T}((x, y), Int[], Int[], T[])
end
const FlexMatrixBuilder{T} = VectorMatrixBuilder{ArrayEntry{T}}
function to_matrix(A::VectorMatrixBuilder)
    sparse(A.Is, A.Js, A.Vs, A.size...)
end
function to_matrix(A::FlexMatrixBuilder)
    _mat = sparse(A.Is, A.Js, A.Vs, A.size..., combine_entries)
    return SparseMatrixCSC(_mat.m, _mat.n, _mat.colptr, _mat.rowval, to_number.(_mat.nzval))
end

function Base.sizehint!(A::VectorMatrixBuilder, n::Int)
    sizehint!(A.Is, n)
    sizehint!(A.Js, n)
    sizehint!(A.Vs, n)
    return A
end
VectorMatrixBuilder{T}(x::Int, y::Int, col_hint::Int) where T =
    sizehint!(VectorMatrixBuilder{T}(x, y), col_hint * x)
Base.@propagate_inbounds function Base.setindex!(A::FlexMatrixBuilder, B::Number, i1::Int, i2::Int; overwrite=true, factor=1)
    # number increment
    push!(A.Is, i1)
    push!(A.Js, i2)
    push!(A.Vs, ArrayEntry(B * factor, overwrite))
    return nothing
end
Base.@propagate_inbounds function Base.setindex!(A::VectorMatrixBuilder, B::Number, i1::Int, i2::Int; overwrite=true, factor=1)
    # number increment
    overwrite && throw(ArgumentError("FastOperatorBuilder does not support overwriting setindex!"))
    push!(A.Is, i1)
    push!(A.Js, i2)
    push!(A.Vs, B * factor)
    return nothing
end
Base.@propagate_inbounds function Base.setindex!(A::AbstractMatrixBuilder, B::AbstractMatrix,
        is1::AbstractArray, is2::AbstractArray; overwrite=true, factor=1)
    for (i1, j1) in enumerate(is1), (i2, j2) in enumerate(is2)
        v = B[i1, i2]
        iszero(v) || (A[j1, j2, overwrite=overwrite, factor=factor] = v)
    end
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

Base.getindex(::AbstractMatrixBuilder, I...) = BuilderView(I)
Base.@propagate_inbounds function Base.setindex!(b::AbstractMatrixBuilder, incr::BuilderIncrementIndex, I...; factor=1)
    @boundscheck incr.axes != I && throw(ArgumentError("Invalid indexing. Only `A[I...] (+=|-=|=) B` is allowed."))
    b[I..., overwrite=false, factor=incr.factor * factor] = incr.B
end

mutable struct UniformMatrixBuilder{T} <: AbstractMatrixBuilder{T}
    sz::Tuple{Int,Int}
    maxcolsize::Int
    colsizes::Vector{Int}
    rowvals::Vector{Int}
    nzvals::Vector{T}
    function UniformMatrixBuilder{T}(szx::Int, szy::Int, col_hint::Int=4) where T
        col_hint = max(col_hint, 1)
        new{T}((szx, szy), col_hint,
            zeros(Int, szy), zeros(Int, col_hint * szy), Array{T}(undef, col_hint * szy))
    end
end

function _grow_to!(b::UniformMatrixBuilder, new_maxcollen)
    new_rowvals = zeros(Int, (new_maxcollen, b.sz[2]))
    new_nzvals = similar(b.nzvals, (new_maxcollen, b.sz[2]))
    copyto!(@view(new_rowvals[1:b.maxcolsize, :]), b.rowvals)
    copyto!(@view(new_nzvals[1:b.maxcolsize, :]), b.nzvals)

    b.rowvals = vec(new_rowvals)
    b.nzvals = vec(new_nzvals)
    b.maxcolsize = new_maxcollen
    return b
end
function Base.sizehint!(b::UniformMatrixBuilder, n::Int)
    new_nrows = n ÷ b.sz[2] + 1
    new_nrows > b.maxcolsize && _grow_to!(b, new_nrows)
    return b
end

Base.@propagate_inbounds function Base.setindex!(b::UniformMatrixBuilder, x::Number, i::Int, j::Int; overwrite=true, factor=1)
    if b.colsizes[j] == b.maxcolsize
        _grow_to!(b, b.maxcolsize * 2)
    end
    col_start = (j - 1) * b.maxcolsize + 1
    col_end = col_start + b.colsizes[j] - 1
    I = col_start
    @inbounds for _ in 1:b.colsizes[j]
        b.rowvals[I] >= i && break
        I += 1
    end
    new_entry = b.rowvals[I] != i
    if b.rowvals[I] != i
        b.colsizes[j] += 1
        @inbounds for i in col_end:-1:I
            b.rowvals[i + 1] = b.rowvals[i]
            b.nzvals[i + 1] = b.nzvals[i]
        end
        b.rowvals[I] = i
    end
    if overwrite || new_entry
        b.nzvals[I] = x * factor
    else
        b.nzvals[I] += x * factor
    end
end

function to_matrix(b::UniformMatrixBuilder)
    mask = b.rowvals .== 0
    rowval = deleteat!(b.rowvals, mask)
    nzval = deleteat!(b.nzvals, mask)
    b.colsizes[1] += 1
    colptr = cumsum(b.colsizes)
    pushfirst!(colptr, 1)
    return SparseMatrixCSC(b.sz[1], b.sz[2], colptr, rowval, nzval)
end

"""
    OperatorBuilder

A helper struct for building custom operators. This struct is used to build operators for a
given system or lattice.
"""
struct OperatorBuilder{SystemT<:System, FieldT<:AbstractField, BuilderT<:AbstractMatrixBuilder}
    sys::SystemT
    mat_builder::BuilderT
    field::FieldT
    auto_hermitian::Bool
end
function OperatorBuilder(sys::SystemT, mat_builder::BuilderT; field::FieldT=NoField(),
        auto_hermitian=false, auto_pbc_field=true) where
    {SystemT<:System, FieldT<:AbstractField, BuilderT<:AbstractMatrixBuilder}
    if auto_pbc_field
        field2 = adapt_field(field, lattice(sys))
    else
        field2 = field
    end
    return OperatorBuilder(sys, mat_builder, field2, auto_hermitian)
end

"""
    OperatorBuilder([T, ]sys, [; field, auto_hermitian, auto_pbc_field])
    OperatorBuilder([T, ]lat, [internal; field, auto_hermitian])

Construct an `OperatorBuilder` for a given system or lattice.

## Arguments
- `T`: The type of the matrix elements. Defaults to `ComplexF64`.
- `sys`: A `System` object representing the system.
- `lat`: The lattice on which the operator is defined.
- `internal`: The basis for the internal degrees of freedom.

## Keyword arguments
- `field`: The gauge field to use for the hopping operators. Defaults to `NoField()`, which
    corresponds to zero magnetic field.
- `auto_hermitian`: Whether to automatically add the hermitian conjugate of the operator.
    Defaults to `false`.
- `auto_pbc_field`: Whether to automatically adapt the field to the periodic boundary
    conditions of the lattice. Defaults to `true`.

## Example
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(5, 5);

julia> builder = OperatorBuilder(l, field=LandauGauge(0.1), auto_hermitian=true)
OperatorBuilder(field=LandauGauge(0.1), auto_hermitian=true)
System: One particle on 25-site SquareLattice in 2D space

julia> hx = Bravais[1, 0]; hy = Bravais[0, 1];

julia> for site in l; builder[site, site + hx] = builder[site, site + hy] = 1; end

julia> H = Hamiltonian(builder)
Hamiltonian(dim=25x25)
System: One particle on 25-site SquareLattice in 2D space
25×25 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 80 stored entries:
⎡⠪⡢⡈⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠢⡈⠠⡢⡈⠢⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠈⠢⡈⠊⡠⡈⠢⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠢⡈⠪⠂⡈⠢⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠢⡈⠪⡢⠈⠢⡀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⠪⡢⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠈⠀⎦

julia> H == tightbinding_hamiltonian(l, field=LandauGauge(0.1))
true
```
"""
function OperatorBuilder(BT::Type{BuilderType}, sys::SystemT; col_hint=nothing, kw...) where
        {BuilderType<:AbstractMatrixBuilder, SystemT<:System}
    oneparticle_len = length(onebodybasis(sys))
    if col_hint === nothing
        builder = BT(oneparticle_len, oneparticle_len)
    else
        builder = BT(oneparticle_len, oneparticle_len, col_hint)
    end
    OperatorBuilder(sys, builder; kw...)
end
OperatorBuilder(T::Type{<:Number}, sys::SystemT; kw...) where {SystemT<:System} =
    OperatorBuilder(FlexMatrixBuilder{T}, sys; kw...)
@accepts_system_t OperatorBuilder

sample(opb::OperatorBuilder) = sample(opb.sys)
const OpBuilderWithInternal = OperatorBuilder{<:System{<:SampleWithInternal}}
const OpBuilderWithoutInternal = OperatorBuilder{<:System{<:SampleWithoutInternal}}

function Base.show(io::IO, mime::MIME"text/plain", opb::OperatorBuilder{Sys,Field,BT}) where {Sys,Field,BT}
    print(io, "OperatorBuilder(",
    opb.field == NoField() ? "" : "field=$(opb.field), ",
    "auto_hermitian=$(opb.auto_hermitian))\nSystem: ")
    show(io, mime, opb.sys)
    if (BT <: VectorMatrixBuilder) && !(BT <: FlexMatrixBuilder)
        print(io, "\nOnly increment/decrement assignments allowed")
    end
end

Base.sizehint!(opb::OperatorBuilder, n::Int) = sizehint!(opb.mat_builder, n)

_internal_one_mat(sample::SampleWithInternal) = internal_one(sample).data
_internal_one_mat(::SampleWithoutInternal) = SMatrix{1,1}(1)
op_to_matrix(sample::Sample, n::Number) = n * _internal_one_mat(sample)
function op_to_matrix(sample::SampleWithInternal, op::DataOperator)
    @boundscheck QuantumOpticsBase.check_samebases(basis(op), internal_basis(sample))
    return op.data
end
function op_to_matrix(sample::Sample, mat::AbstractMatrix)
    @boundscheck all(==(internal_length(sample)), size(mat)) ||
        throw(ArgumentError("matrix size does not match on-site dims"))
    return mat
end
op_to_matrix(::SampleWithoutInternal, ::DataOperator) =
    throw(ArgumentError("cannot define on-site operator for a system without internal degrees of freedom"))
op_to_matrix(::Sample, op) =
    throw(ArgumentError("cannot interpret $(typeof(op)) as on-site operator"))

const UnresolvedSite = Union{AbstractSite, ResolvedSite, Int}
function expand_bond(l::AbstractLattice, site1::UnresolvedSite, site2::UnresolvedSite, field::AbstractField)
    s1 = resolve_site(l, site1)
    s2 = resolve_site(l, site2)
    s1 === nothing && return nothing
    s2 === nothing && return nothing
    return expand_bond(l, s1, s2, field)
end

@inline function expand_bond(::AbstractLattice, s1::ResolvedSite, s2::ResolvedSite, field::AbstractField)
    field_fact = exp(-2π * im * line_integral(field, s1.old_site.coords, s2.old_site.coords))
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

function _construct_manybody_maybe(sys::ManyBodySystem, op::AbstractOperator)
    bas = basis(op)
    check_samebases(bas, onebodybasis(sys))
    return manybodyoperator(ManyBodyBasis(bas, occupations(sys)), op)
end
_construct_manybody_maybe(::OneParticleBasisSystem, op::AbstractOperator) = op
function to_operator(opb::OperatorBuilder, additional_term=nothing)
    op = Operator(onebodybasis(opb.sys), to_matrix(opb.mat_builder))
    if additional_term === nothing
        return _construct_manybody_maybe(opb.sys, op)
    else
        return _construct_manybody_maybe(opb.sys, op + additional_term)
    end
end
QuantumOpticsBase.Operator(opb::OperatorBuilder) = to_operator(opb)
Hamiltonian(opb::OperatorBuilder) = Hamiltonian(opb.sys, Operator(opb))
