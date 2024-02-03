import QuantumOpticsBase: DataOperator

function add_diagonal_term!(builder, op, diag)
    for i in eachindex(diag)
        increment!(builder, op, i, i, factor=diag[CartesianIndex(i)])
    end
end
function add_diagonal_term!(builder, op)
    for i in eachindex(diag)
        increment!(builder, op, i, i)
    end
end

function add_hopping_term!(builder, l::AbstractLattice, op, bond::SiteOffset,
        field::AbstractField)
    dims(bond) > dims(l) && error("Incompatible dims")
    for site1 in l
        site2 = site1 + bond
        site2 === NoSite() && continue
        add_hopping_term!(builder, l, op, site1 => site2, field)
    end
end
function add_hopping_term!(builder, l::AbstractLattice, op, bond::SingleBond,
        field::AbstractField)
    prop = expand_bond(l, bond, field)
    prop === nothing && return
    i, j, total_factor = prop
    !isfinite(total_factor) && error("got NaN or Inf when finding the phase factor")
    increment!(builder, op, i, j, factor=total_factor)
    increment!(builder, op', j, i, factor=total_factor')
end

const AbstractBond = Union{SiteOffset, SingleBond}
function add_term!(builder::AbstractMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:Tuple}; kw...)
    # Expand 'lattice' part
    opdata, lparts = arg
    for lpart in lparts
        add_term!(builder, sample, opdata => lpart; kw...)
    end
end
function add_term!(builder::AbstractMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractBond};
        field=NoField())
    # Hopping term
    opdata, bond = arg
    add_hopping_term!(builder, sample.latt, opdata, bond, field)
end
function add_term!(builder::AbstractMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:LatticeValue}; kw...)
    # Diagonal operator
    opdata, lv = arg
    add_diagonal_term!(builder, opdata, lv.values)
end
function add_term!(builder::AbstractMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:Nothing}; kw...)
    # Diagonal operator, lattice part ==1
    add_diagonal_term!(builder, arg[1])
end
function add_term!(builder::AbstractMatrixBuilder, ::Sample, arg::DataOperator; kw...)
    # Arbitrary sparse operator
    increment!(builder, arg.data)
end

function arg_to_pair(sample::Sample, arg::DataOperator)
    if samebases(basis(arg), basis(sample))
        return sparse(arg)
    elseif hasinternal(sample) && samebases(basis(arg), internal_basis(sample))
        return arg.data => nothing
    elseif samebases(basis(arg), LatticeBasis(lattice(sample)))
        if hasinternal(sample)
            return sparse(arg)
        else
            return internal_one(sample) ⊗ sparse(arg)
        end
    else
        error("Invalid Operator argument: basis does not match neither lattice nor internal phase space")
    end
end

function arg_to_pair(sample::Sample, arg::AbstractMatrix)
    if all(==(internal_length(sample)), size(arg))
        return arg => nothing
    else
        error("Invalid Matrix argument: size does not match on-site dimension count")
    end
end
arg_to_pair(sample::Sample, n::Number) = op_to_matrix(sample, n) => nothing
arg_to_pair(sample::Sample, arg) = _internal_one_mat(sample) => arg

struct Nearest{N}
    is::NTuple{N, Int}
    Nearest(is::Vararg{Int, N}) where N = new{N}(is)
end
get_bonds(l::AbstractLattice, n::Nearest) = tuple(((default_bonds(l, i) for i in n.is)...)...)
function arg_to_pair(sample::Sample, arg::Pair)
    op, lpart = arg
    if lpart isa LatticeValue
        check_samesites(lpart, sample)
        new_lpart = lpart
    elseif lpart isa AbstractBond || lpart isa Tuple{Vararg{<:AbstractBond}}
        new_lpart = lpart
    elseif lpart isa Nearest
        new_lpart = get_bonds(lattice(sample), lpart)
    else
        throw(ArgumentError("unsupported on-lattice operator type $(typeof(lpart))"))
    end
    return op_to_matrix(sample, op) => new_lpart
end

function build_operator(T::Type, sys::System, args...; field=NoField())
    sample = sys.sample
    builder = SparseMatrixBuilder{T}(length(sample), length(sample))
    for arg in args
        add_term!(builder, sample, arg_to_pair(sample, arg);
            field=field)
    end
    op = Operator(basis(sample), to_matrix(builder))
    return _build_manybody_maybe(sys, op)
end
@accepts_system_t build_operator
build_hamiltonian(T::Type, sys::System, args...; kw...) =
    Hamiltonian(sys, build_operator(T, sys, args...; kw...))
@accepts_system_t build_hamiltonian

tightbinding_hamiltonian(T::Type, sys::System; t1=1, t2=0, t3=0, field=NoField()) =
    build_hamiltonian(T, sys, t1 => Nearest(1), t2 => Nearest(2), t3 => Nearest(3), field=field)
@accepts_system_t tightbinding_hamiltonian
