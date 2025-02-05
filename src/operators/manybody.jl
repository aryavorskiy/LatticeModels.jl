import QuantumOpticsBase: Basis, SparseOperator
import QuantumOpticsBase: allocate_buffer, state_index, state_transition!

_count_onsite(occ::AbstractVector, i, N) = sum(@view occ[(i - 1) * N + 1: i * N])
function _count_onsite(occ::FermionBitstring, i, N)
    mask1 = (one(occ.bits) << i * N) - 1
    mask2 = (one(occ.bits) << (i - 1) * N) - 1
    return count_ones(occ.bits & mask1 & ~mask2)
end

function _2p_interaction_collect(f::Function, T::Type{<:Number}, lat::AbstractLattice)
    M = zeros(T, length(lat), length(lat))
    for i in eachindex(lat)
        site1 = lat[i]
        for j in 1:i - 1
            site2 = translate_to_nearest(lat, site1, lat[j])
            M[i, j] = f(site1, site2)
            M[j, i] = M[i, j]
        end
        M[i, i] = f(site1, site1)
    end
    return M
end

function _2p_interaction_diags(M::AbstractMatrix, occups, N::Int)
    diags = eltype(M)[]
    for occ in occups
        int_energy = 0.
        for i in 1:size(M, 1)
            occi = _count_onsite(occ, i, N)
            occi == 0 && continue
            for j in 1:i - 1
                occj = _count_onsite(occ, j, N)
                occj == 0 && continue
                int_energy += M[i, j] * occi * occj
            end
            int_energy += M[i, i] * (occi * (occi - 1) ÷ 2)
        end
        push!(diags, int_energy)
    end
    return diags
end

"""
    interaction(f, [T, ]sys)

Create an two-site interaction operator for a given `NParticles` system. The function `f` takes two
arguments, which are the two sites, and returns the interaction energy.
"""
function interaction(f::Function, T::Type{<:Number}, sys::NParticles; occupations_type = nothing)
    l = lattice(sys)
    M = _2p_interaction_collect(f, T, l)
    occups = occupations(sys, occupations_type)
    diags = _2p_interaction_diags(M, occups, internal_length(sys))
    diagonaloperator(ManyBodyBasis(onebodybasis(sys), occups), diags)
end

"""
    interaction(f, [T, ]sys, K)

Create an `2K`-site interaction operator for a given `NParticles` system. The function `f`
takes two `K`-tuples of integer numbers, which are site indices for creation and annihilation
operators, and returns the interaction energy.

If the system `sys` has internal degrees of freedom, the function `f` should take four `K`-tuples:
first two are site & internal indices for creation operators, and the last two are the same for
annihilation operators.
"""
function interaction(f::Function, T::Type{<:Number}, sys::NParticles, ::Val{K}) where K
    # This function uses undocumented QuantumOpticsBase API. Be careful!
    @assert K ≤ sys.nparticles
    l = lattice(sys)
    N = internal_length(sys)
    le = length(l) * N
    c_indices = CartesianIndex{K}():CartesianIndex{K}(fill(le, K)...)
    l_inds(inds) = inds .|> (i -> l[(i - 1) ÷ N + 1])
    i_inds(inds) = inds .|> (i -> (i - 1) % N + 1)

    is = Int[]
    js = Int[]
    vs = T[]

    bas = basis(sys)
    buffer = allocate_buffer(bas)
    for a_inds_c in c_indices
        a_inds = Tuple(a_inds_c)
        issorted(a_inds) || continue
        for at_inds_c in c_indices
            at_inds = Tuple(at_inds_c)
            issorted(at_inds) || continue
            for (m, occ) in enumerate(bas.occupations)
                C = state_transition!(buffer, occ, at_inds, a_inds)
                C === nothing && continue
                value = if hasinternal(sys)
                    f(l_inds(at_inds), i_inds(at_inds), l_inds(a_inds), i_inds(a_inds))::Number
                else
                    f(l_inds(at_inds), l_inds(a_inds))::Number
                end
                value == 0 && continue
                n = state_index(bas.occupations, buffer)
                push!(is, m)
                push!(js, n)
                push!(vs, C * value)
            end
        end
    end
    occ_len = length(bas.occupations)
    return SparseOperator(bas, sparse(is, js, vs, occ_len, occ_len))
end
interaction(f::Function, T::Type{<:Number}, sys::NParticles, K::Int; kw...) =
    interaction(f, T, sys, Val(K); kw...)
interaction(f::Function, args...; kw...) = interaction(f, ComplexF64, args...; kw...)
interaction(f::Function, ::Type, ::Type, args...; kw...) = throw(MethodError(interaction, (f, args...)))
