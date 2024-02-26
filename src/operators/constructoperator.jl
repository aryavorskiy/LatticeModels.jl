import QuantumOpticsBase: DataOperator

function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:SingleBond})
    # Single bond
    op, bond = arg
    builder[bond[1], bond[2]] += op
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:AbstractBonds})
    # Hopping term
    op, bonds = arg
    for (s1, s2) in bonds
        builder[s1, s2] += op
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:DirectedBonds})
    # Hopping term
    op, bonds = arg
    lat = lattice(builder)
    for i in 1:length(lat)
        site = lat[i]
        s1 = ResolvedSite(site, i)
        for site2 in destinations(bonds, site)
            s2 = resolve_site(lat, site2)
            s2 !== nothing && (builder[s1, s2] += op)
        end
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:LatticeValue})
    # Diagonal operator
    op, lv = arg
    N = internal_length(builder)
    for i in 1:length(lattice(builder)) # This avoids finding indices of sites
        is = (i - 1) * N + 1:i * N
        builder.mat_builder[is, is, factor=lv.values[i], overwrite=false] = op
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:Number})
    # Diagonal operator, lattice part ==1
    op, n = arg
    N = internal_length(builder)
    for i in 1:length(lattice(builder)) # This avoids finding indices of sites
        is = (i - 1) * N + 1:i * N
        builder.mat_builder[is, is, factor=n, overwrite=false] = op
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
        throw(ArgumentError("operator basis does not match neither lattice nor internal phase space"))
    end
end

arg_to_pair(sample::Sample, arg::AbstractMatrix) = op_to_matrix(sample, arg) => 1
arg_to_pair(sample::Sample, arg) = _internal_one_mat(sample) => arg

function arg_to_pair(sample::Sample, arg::Pair)
    op, lpart = arg
    if lpart isa LatticeValue
        check_samesites(lpart, sample)
        new_lpart = lpart
    elseif lpart isa AbstractBonds
        new_lpart = adapt_bonds(lpart, lattice(sample))
    elseif lpart isa Union{Number, SingleBond}
        new_lpart = lpart
    else
        throw(ArgumentError("cannot interpret $(typeof(lpart)) as on-lattice operator"))
    end
    return op_to_matrix(sample, op) => new_lpart
end

function construct_operator(T::Type, sys::System, args...; field=NoField())
    sample = sys.sample
    builder = FastOperatorBuilder(T, sys; field=field, auto_hermitian=true)
    sizehint!(builder.mat_builder, length(sample) * length(args))
    for arg in args
        pair = arg_to_pair(sample, arg)
        iszero(first(pair)) || add_term!(builder, pair)
    end
    Operator(builder)
end
@accepts_system_t construct_operator

construct_hamiltonian(T::Type, sys::System, args...; kw...) =
    Hamiltonian(sys, construct_operator(T, sys, args...; kw...))
@accepts_system_t construct_hamiltonian

tightbinding_hamiltonian(T::Type, sys::System, args...; t1=1, t2=0, t3=0, field=NoField()) =
    construct_hamiltonian(T, sys, args...,
        t1 => NearestNeighbor(1),
        t2 => NearestNeighbor(2),
        t3 => NearestNeighbor(3), field=field)
@accepts_system_t tightbinding_hamiltonian
