import QuantumOpticsBase: DataOperator

struct SparseMatrixBuilder{T} <: AbstractMatrix{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    SparseMatrixBuilder{T}(sz) where T = new{T}(sz, [], [], [])
end
SparseMatrixBuilder{T}(sz...) where T = SparseMatrixBuilder{T}(sz)
SparseMatrixBuilder{T}(size_x::Int) where T = SparseMatrixBuilder{T}(size_x, size_x)
Base.size(smb::SparseMatrixBuilder) = smb.size
to_matrix(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int; factor=1)
    @boundscheck @assert 1 ≤ i1 ≤ builder.size[1]
    @boundscheck @assert 1 ≤ i2 ≤ builder.size[2]
    push!(builder.Is, i1)
    push!(builder.Js, i2)
    push!(builder.Vs, rhs * factor)
end

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::AbstractMatrix, i1::Int, i2::Int; factor=1)
    N = size(rhs)[1]
    for i in 1:N, j in 1:N
        v = rhs[i, j]
        iszero(v) || increment!(builder, v, i + N * (i1 - 1), j + N * (i2 - 1); factor=factor)
    end
end
Base.@propagate_inbounds increment!(builder::SparseMatrixBuilder, rhs::DataOperator, i1, i2) =
    increment!(builder, rhs.data, i1, i2)

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::SparseMatrixCSC; factor=1)
    @boundscheck @assert size(rhs) == size(builder)
    nis, njs, nvs = findnz(rhs)
    append!(builder.Is, nis)
    append!(builder.Js, njs)
    append!(builder.Vs, nvs * factor)
    nothing
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

struct OperatorBuilder{SystemT, FieldT, T}
    sys::SystemT
    field::FieldT
    internal_length::Int
    builder::SparseMatrixBuilder{T}
    auto_hermitian::Bool
    OperatorBuilder{T}(sys::SystemT, field::FieldT=NoField(); auto_hermitian=false) where {SystemT<:System, FieldT<:AbstractField, T} =
        new{SystemT, FieldT, T}(sys, field, internal_length(sys), SparseMatrixBuilder{T}(length(onebodybasis(sys))), auto_hermitian)
end
@accepts_lattice OperatorBuilder
lattice(opb::OperatorBuilder) = lattice(opb.sys)
OperatorBuilder(args...; kw...) = OperatorBuilder{ComplexF64}(args...; kw...)
const OpBuilderWithInternal = OperatorBuilder{<:System{<:SampleWithInternal}}
const OpBuilderWithoutInternal = OperatorBuilder{<:System{<:SampleWithoutInternal}}

check_arg(::OpBuilderWithoutInternal, n::Number) = n
check_arg(opb::OperatorBuilder, n::Number) = n * internal_one(opb.sys)
Base.@propagate_inbounds function check_arg(opb::OpBuilderWithInternal, op::DataOperator)
    @boundscheck QuantumOpticsBase.check_samebases(basis(op), internal_basis(opb.sys))
    return op.data
end
Base.@propagate_inbounds function check_arg(opb::OpBuilderWithInternal, mat::AbstractMatrix)
    @boundscheck @assert all(==(internal_length(opb.sys)), size(mat))
    return mat
end

Base.@propagate_inbounds function increment!(opbuilder::OperatorBuilder, rhs, lp1, lp2)
    new_rhs = check_arg(opbuilder, rhs)
    l = lattice(opbuilder)
    site1 = get_site(l, lp1)
    site2 = get_site(l, lp2)
    ifact, new_site1 = shift_site(opbuilder.sys, site1)
    jfact, new_site2 = shift_site(opbuilder.sys, site2)
    field_fact = exp(-2π * im * line_integral(opbuilder.field, site1.coords, site2.coords))
    l = lattice(opbuilder)
    i1 = site_index(l, new_site1)
    i2 = site_index(l, new_site2)
    i1 === nothing && return
    i2 === nothing && return
    total_factor = jfact * field_fact * ifact'
    increment!(opbuilder.builder, new_rhs, i1, i2, factor = total_factor)
    if opbuilder.auto_hermitian && i1 != i2
        increment!(opbuilder.builder, new_rhs', i2, i1, factor = total_factor')
    end
end

function to_operator(opb::OperatorBuilder; warning=false)
    op = Operator(onebodybasis(opb.sys), to_matrix(opb.builder))
    if !ishermitian(op) && !warning && !opb.auto_hermitian
        @warn "The resulting operator is not hermitian. Set `warning=false` to hide this message or add `auto_hermitian=true` to the `OperatorBuilder` constructor."
    end
    return manybodyoperator(opb.sys, op)
end

struct Hamiltonian{SystemT, BasisT, T} <: DataOperator{BasisT, BasisT}
    sys::SystemT
    basis_l::BasisT
    basis_r::BasisT
    data::T
end
function Hamiltonian(sys::System, op::Operator)
    return Hamiltonian(sys, basis(op), basis(op), op.data)
end
QuantumOpticsBase.Operator(ham::Hamiltonian) = Operator(ham.basis_l, ham.data)
system(::DataOperator) = nothing
system(ham::Hamiltonian) = ham.sys

function tightbinding_hamiltonian(sys::System; t1=1, t2=0, t3=0,
    field=NoField(), boundaries=BoundaryConditions())
    sample = sys.sample
    l = lattice(sys)
    builder = SparseMatrixBuilder{ComplexF64}(length(sys), length(sys))
    internal_eye = internal_one(sample)
    for bond in default_bonds(l)
        add_hoppings!(builder, nothing, l, t1 * internal_eye, bond, field, boundaries)
    end
    if t2 != 0
        for bond in default_nnbonds(l)
            add_hoppings!(builder, nothing, l, t2 * internal_eye, bond, field, boundaries)
        end
    end
    if t3 != 0
        for bond in default_nnnbonds(l)
            add_hoppings!(builder, nothing, l, t3 * internal_eye, bond, field, boundaries)
        end
    end
    return Hamiltonian(sys, manybodyoperator(sys, to_matrix(builder)))
end
@accepts_lattice tightbinding_hamiltonian

const AbstractSiteOffset = Union{SiteOffset, SingleBond}
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractSiteOffset};
        field=NoField())
    # Hopping operator
    opdata, bond = arg
    add_hoppings!(builder, sample.adjacency_matrix, sample.latt, opdata, bond, field, sample.boundaries)
end
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:LatticeValue}; kw...)
    # Diagonal operator
    opdata, lv = arg
    add_diagonal!(builder, opdata, lv.values)
end
function build_operator!(builder::SparseMatrixBuilder, ::Sample, arg::DataOperator; kw...)
    # Arbitrary sparse operator
    increment!(builder, arg.data)
end

function preprocess_argument(sample::Sample, arg::DataOperator)
    if samebases(basis(arg), onebodybasis(sample))
        sparse(arg)
    elseif samebases(basis(arg), sample.internal)
        sparse(arg) ⊗ one(LatticeBasis(sample.latt))
    elseif samebases(basis(arg), LatticeBasis(sample.latt))
        internal_one(sample) ⊗ sparse(arg)
    else
        error("Invalid Operator argument: basis does not match neither lattice nor internal phase space")
    end
end

function preprocess_argument(sample::Sample, arg::AbstractMatrix)
    bas = sample.internal
    if all(==(length(bas)), size(arg))
        preprocess_argument(sample, SparseOperator(bas, arg))
    else
        error("Invalid Matrix argument: size does not match on-site dimension count")
    end
end

preprocess_argument(sample::Sample, arg) =
    preprocess_argument(sample, internal_one(sample) => arg)

preprocess_argument(sample::Sample, n::Number) = preprocess_argument(sample, n * internal_one(sample))

function preprocess_argument(sample::Sample, arg::Pair)
    op, on_lattice = arg
    if op isa Operator
        check_samebases(basis(op), sample.internal)
        opdata = sparse(op.data)
    elseif op isa AbstractMatrix
        opdata = sparse(op)
    elseif op isa Number
        opdata = op
    else
        error("Invalid Pair argument: unsupported on-site operator type")
    end
    if on_lattice isa LatticeValue
        check_lattice_match(on_lattice, sample)
        return opdata => on_lattice
    elseif on_lattice isa Number
        return preprocess_argument(sample, opdata * on_lattice)
    elseif on_lattice isa AbstractSiteOffset
        return opdata => on_lattice
    else
        error("Invalid Pair argument: unsupported on-lattice operator type")
    end
end

function build_hamiltonian(sys::System, args...;
        field=NoField())
    sample = sys.sample
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    for arg in args
        build_operator!(builder, sample, preprocess_argument(sample, arg);
            field=field)
    end
    op = Operator(onebodybasis(sample), to_matrix(builder))
    return Hamiltonian(sys, manybodyoperator(sys, op))
end
@accepts_lattice build_hamiltonian

hoppings(adj, l::Lattice, bs::SiteOffset...; kw...) = Operator(build_hamiltonian(Sample(adj, l), bs...; kw...))
hoppings(l::Lattice, bs::SiteOffset...; kw...) = Operator(build_hamiltonian(Sample(l), bs...; kw...))
