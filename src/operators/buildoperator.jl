import QuantumOpticsBase: DataOperator

function add_hopping_term!(builder::OperatorBuilder, op, bonds::SiteOffset)

end

const AbstractBond = Union{SiteOffset, SingleBond}
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:Tuple})
    # Expand 'lattice' part
    opdata, lparts = arg
    for lpart in lparts
        add_term!(builder, opdata => lpart)
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:SingleBond})
    # Hopping term
    op, bond = arg
    @increment builder[bond[1], bond[2]] = op
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:SiteOffset})
    # Hopping term
    op, bonds = arg
    dims(bonds) > dims(builder) && error("Incompatible dims")
    for site1 in lattice(builder)
        @increment builder[site1, site1 + bonds] += op
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:LatticeValue})
    # Diagonal operator
    op, lv = arg
    N = internal_length(builder)
    for i in 1:length(lattice(builder)) # This avoids finding indices of sites
        increment!(builder.mat_builder, op, (i-1)*N+1:i*N, (i-1)*N+1:i*N, factor=lv.values[i])
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:Number})
    # Diagonal operator, lattice part ==1
    op, n = arg
    N = internal_length(builder)
    for i in 1:length(lattice(builder)) # This avoids finding indices of sites
        increment!(builder.mat_builder, op, (i-1)*N+1:i*N, (i-1)*N+1:i*N, factor=n)
    end
end
function add_term!(builder::OperatorBuilder, arg::DataOperator)
    # Arbitrary sparse operator
    increment!(builder.mat_builder, arg.data)
end

function arg_to_pair(sample::Sample, arg::DataOperator)
    if samebases(basis(arg), basis(sample))
        return sparse(arg)
    elseif hasinternal(sample) && samebases(basis(arg), internal_basis(sample))
        return arg.data => 1
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
        return arg => 1
    else
        error("Invalid Matrix argument: size does not match on-site dimension count")
    end
end
arg_to_pair(sample::Sample, n::Number) = op_to_matrix(sample, n) => 1
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
    elseif lpart isa AbstractBond || lpart isa Tuple{Vararg{AbstractBond}}
        new_lpart = lpart
    elseif lpart isa Nearest
        new_lpart = get_bonds(lattice(sample), lpart)
    elseif lpart isa Number
        new_lpart = lpart
    else
        throw(ArgumentError("unsupported on-lattice operator type $(typeof(lpart))"))
    end
    return op_to_matrix(sample, op) => new_lpart
end

function build_operator(T::Type, sys::System, args...; field=NoField())
    sample = sys.sample
    builder = FastSparseOperatorBuilder(T, sys; field=field, auto_hermitian=true)
    for arg in args
        add_term!(builder, arg_to_pair(sample, arg))
    end
    Operator(builder)
end
@accepts_system_t build_operator
build_hamiltonian(T::Type, sys::System, args...; kw...) =
    Hamiltonian(sys, build_operator(T, sys, args...; kw...))
@accepts_system_t build_hamiltonian

tightbinding_hamiltonian(T::Type, sys::System; t1=1, t2=0, t3=0, field=NoField()) =
    build_hamiltonian(T, sys, t1 => Nearest(1), t2 => Nearest(2), t3 => Nearest(3), field=field)
@accepts_system_t tightbinding_hamiltonian
