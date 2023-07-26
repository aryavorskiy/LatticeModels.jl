struct SparseMatrixBuilder{T} <: AbstractMatrix{T}
    size::Tuple{Int,Int}
    Is::Vector{Int}
    Js::Vector{Int}
    Vs::Vector{T}
    SparseMatrixBuilder{T}(sz) where T = new{T}(sz, [], [], [])
    SparseMatrixBuilder{T}(sz...) where T = SparseMatrixBuilder{T}(sz)
end
Base.size(smb::SparseMatrixBuilder) = smb.size
to_matrix(builder::SparseMatrixBuilder) =
    sparse(builder.Is, builder.Js, builder.Vs, builder.size...)

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::Number, i1::Int, i2::Int, factor=1)
    @boundscheck @assert 1 ≤ i1 ≤ builder.size[1]
    @boundscheck @assert 1 ≤ i2 ≤ builder.size[2]
    push!(builder.Is, i1)
    push!(builder.Js, i2)
    push!(builder.Vs, rhs * factor)
end

Base.@propagate_inbounds function increment!(builder::SparseMatrixBuilder, rhs::AbstractMatrix, i1::Int, i2::Int, factor=1)
    N = size(rhs)[1]
    for i in 1:N, j in 1:N
        v = rhs[i, j]
        iszero(v) || increment!(builder, v, i + N * (i1 - 1), j + N * (i2 - 1), factor)
    end
end

function increment!(builder::SparseMatrixBuilder, rhs::SparseMatrixCSC)
    @assert size(rhs) == size(builder)
    nis, njs, nvs = findnz(rhs)
    append!(builder.Is, nis)
    append!(builder.Js, njs)
    append!(builder.Vs, nvs)
    nothing
end

function tightbinding_hamiltonian(sample::Sample; tn=1, tnn=0, tnnn=0,
    field=NoField(), boundaries=BoundaryConditions())
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    internal_eye = one(internal).data
    for bond in default_bonds(l)
        add_hoppings!(builder, nothing, l, tn * internal_eye, bond, field, boundaries)
    end
    if tnn != 0
        for bond in default_nnbonds(l)
            add_hoppings!(builder, nothing, l, tnn * internal_eye, bond, field, boundaries)
        end
    end
    if tnnn != 0
        for bond in default_nnnbonds(l)
            add_hoppings!(builder, nothing, l, tnnn * internal_eye, bond, field, boundaries)
        end
    end
    return manybodyoperator(sample, to_matrix(builder))
end
@accepts_lattice tightbinding_hamiltonian

const AbstractBonds = Union{Bonds, SingleBond}
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractBonds};
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
function build_operator!(builder::SparseMatrixBuilder, ::Sample, arg::Operator; kw...)
    # Arbitrary sparse operator
    increment!(builder, arg.data)
end

function preprocess_argument(sample::Sample, arg::Operator)
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
    elseif on_lattice isa AbstractBonds
        return opdata => on_lattice
    else
        error("Invalid Pair argument: unsupported on-lattice operator type")
    end

end

function build_hamiltonian(sample::Sample, args...;
        field=NoField())
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    for arg in args
        build_operator!(builder, sample, preprocess_argument(sample, arg);
            field=field)
    end
    op = Operator(onebodybasis(sample), to_matrix(builder))
    return manybodyoperator(sample, op)
end
@accepts_lattice build_hamiltonian

hoppings(adj, l::Lattice, bs::Bonds...; kw...) = build_hamiltonian(Sample(adj, l), bs...; kw...)
hoppings(l::Lattice, bs::Bonds...; kw...) = build_hamiltonian(Sample(l), bs...; kw...)
