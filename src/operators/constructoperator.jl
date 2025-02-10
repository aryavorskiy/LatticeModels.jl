import QuantumOpticsBase: DataOperator
using LinearAlgebra

function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:SingleBond})
    # Single bond
    op, bond = arg
    builder[bond[1], bond[2]] += op
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:AbstractSite})
    # Single bond
    op, site = arg
    builder[site, site] += op
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:AbstractBonds})
    # Hopping term
    op, bonds = arg
    for (s1, s2) in bonds
        builder[s1, s2] += op
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:LatticeValue})
    # Onsite term
    op, lv = arg
    N = internal_length(builder)
    Is = 1:length(lattice(builder))
    for i in Is    # This avoids finding indices of sites
        is = (i - 1) * N + 1:i * N
        builder.mat_builder[is, is, factor=lv.values[i], overwrite=false] = op
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:BravaisSiteMapping})
    # Bravais site mapping
    l = lattice(builder)
    op, bsm = arg
    for translation in bsm.translations
        add_term!(builder, op => adapt_bonds(translation, l))
    end
end
function add_term!(builder::OperatorBuilder, arg::Pair{<:Any, <:Number})
    # Onsite term, lattice part ==1
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

function add_terms!(opb::OperatorBuilder, arg, args...)
    pair = arg_to_pair(sample(opb), arg)
    (pair isa Pair && iszero(first(pair))) || add_term!(opb, pair)
    add_terms!(opb, args...)
end
add_terms!(::OperatorBuilder) = nothing

function arg_to_pair(sample::Sample, arg::DataOperator)
    if samebases(basis(arg), basis(sample))
        return sparse(arg)
    elseif hasinternal(sample) && samebases(basis(arg), internal_basis(sample))
        return arg.data => 1
    elseif samebases(basis(arg), LatticeBasis(lattice(sample)))
        if hasinternal(sample)
            return internal_one(sample) ⊗ sparse(arg)
        else
            return sparse(arg)
        end
    else
        throw(ArgumentError("operator basis does not match neither lattice nor internal phase space"))
    end
end

arg_to_pair(sample::Sample, arg::AbstractMatrix) = op_to_matrix(sample, arg) => 1
arg_to_pair(sample::Sample, arg) = _internal_one_mat(sample) => arg

function arg_to_pair(sample::Sample, arg::Pair{<:Union{Number,AbstractMatrix,DataOperator}})
    op, onlat = arg
    if onlat isa LatticeValue
        check_samesites(onlat, sample)
        new_onlat = onlat
    elseif onlat isa AbstractBonds
        new_onlat = adapt_bonds(onlat, lattice(sample))
    elseif onlat isa Union{Number, SingleBond, AbstractSite}
        new_onlat = onlat
    else
        throw(ArgumentError("cannot interpret $(typeof(onlat)) as on-lattice operator"))
    end
    return op_to_matrix(sample, op) => new_onlat
end

_colsize_est(p::Pair) = _colsize_est(p[1]) * _colsize_est(p[2])
# onlat
_colsize_est(::Number) = 1
_colsize_est(::LatticeValue) = 1
_colsize_est(::AbstractTranslation) = 2
_colsize_est(btr::BravaisTranslation) = 1 + allequal(btr.site_indices)
_colsize_est(bsm::BravaisSiteMapping) = sum(_colsize_est, bsm.translations)
_colsize_est(::AbstractBonds) = 4
_colsize_est(::AbstractSite) = 0
_colsize_est(::SingleBond) = 1
# op
_colsize_est(mat::SparseMatrixCSC) = maximum(diff(mat.colptr))
_colsize_est(mat::AbstractMatrix) = maximum(count(!=(0), mat, dims=1))
_colsize_est(::Diagonal) = 1
_colsize_est(op::DataOperator) = _colsize_est(op.data)

"""
    construct_operator([T, ]sys, terms...[; field, auto_pbc_field])
    construct_operator([T, ]lat[, internal, terms...; field])

Construct an operator for the given system.

Each of the `terms` describes a term of the Hamiltonian. The term can be given in several ways:
- A `DataOperator` on the lattice, internal or composite basis (will be matched automatically).
- A `Pair` of an "internal" and an "on-lattice" part (e.g. `int_p => lat_p`):
    - The "internal" part can be a `DataOperator`, a matrix or a number.
    - The "on-lattice" part can be a `LatticeValue` (represents a diagonal term), a site
        (represents a local on-site potential), a bond (represents a hopping term) or a `site1 => site2`
        pair (represents a single hopping).
    - Identity "internal" or "on-lattice" parts can be omitted.

See documentation for more details.

## Arguments
- `T`: The element type of the Hamiltonian. Default is `ComplexF64`.
- `sys`: The `System` for which the Hamiltonian is constructed.
- `lat`: The lattice for which the Hamiltonian is constructed.
- `internal`: The basis for the internal degrees of freedom.

## Keyword Arguments
- `field`: The gauge field to use for the bond operators. Default is `NoField()`.
- `auto_pbc_field`: Whether to automatically adapt the field to the periodic boundary
    conditions of the lattice. Defaults to `true`.
"""
function construct_operator(BT::Type{<:AbstractMatrixBuilder{T}}, sys::System, args...; kw...) where T
    opb = OperatorBuilder(BT, sys; auto_hermitian=true, kw...)
    add_terms!(opb, args...)
    return Operator(opb)
end
function construct_operator(T::Type{<:Number}, sys::System, args...; kw...)
    arg_pairs = map(arg -> arg_to_pair(sys.sample, arg), args)
    colsize_hint = sum(_colsize_est, arg_pairs)
    if colsize_hint < 1000
        construct_operator(UniformMatrixBuilder{T}, sys, args...;
            size_hint=colsize_hint, kw...)
    else
        construct_operator(SimpleMatrixBuilder{T}, sys, args...;
            size_hint=colsize_hint, kw...)
    end
end
@accepts_system_t construct_operator

"""
    construct_hamiltonian([T, ]sys, terms...[; field, auto_pbc_field])
    construct_hamiltonian([T, ]lat[, internal, terms...; field, auto_pbc_field])

Construct a Hamiltonian for the given system. Does the same as `construct_operator`, but wraps
the result in a `Hamiltonian` type.
"""
construct_hamiltonian(T::Type, sys::System, args...; kw...) =
    Hamiltonian(sys, construct_operator(T, sys, args...; kw...))
@accepts_system_t construct_hamiltonian

"""
    tightbinding_hamiltonian([T, ]sys[, args...; t1=1, t2=0, t3=0, field, auto_pbc_field])
    tightbinding_hamiltonian([T, ]lat[, internal, args...; t1=1, t2=0, t3=0, field, auto_pbc_field])

Construct a tight-binding Hamiltonian for the given system.

## Arguments
- `T`: The element type of the Hamiltonian. Default is `ComplexF64`.
- `sys`: The `System` for which the Hamiltonian is constructed.
- `lat`: The lattice for which the Hamiltonian is constructed.
- `internal`: The basis for the internal degrees of freedom.
All other arguments are interpreted as terms of the Hamiltonian and passed to `construct_hamiltonian`.

## Keyword Arguments
- `t1`, `t2`, `t3`: The hopping amplitudes for the nearest, next-nearest, and next-next-nearest
    neighbors, respectively.
- `field`: The gauge field to use for the bond operators. Default is `NoField()`.
- `auto_pbc_field`: Whether to automatically adapt the field to the periodic boundary
    conditions of the lattice. Defaults to `true`.
"""
tightbinding_hamiltonian(args...; t1=1, t2=0, t3=0, kw...) =
    construct_hamiltonian(args...,
        t1 => NearestNeighbor(1),
        t2 => NearestNeighbor(2),
        t3 => NearestNeighbor(3); kw...)
