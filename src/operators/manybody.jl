import QuantumOpticsBase: Basis, SparseOperator

"""
    interaction(f, [T, ]sys)

Create an two-site interaction operator for a given `NParticles` system. The function `f` takes two
arguments, which are the two sites, and returns the interaction energy.
"""
function interaction(f::Function, T::Type{<:Number}, sys::NParticles)
    l = lattice(sys)
    occups = occupations(sys)
    diags = T[]
    N = internal_length(sys)
    for occ in occups
        int_energy = 0.
        for i in eachindex(l)
            occi = sum(@view occ[(i - 1) * N + 1: i * N])
            occi == 0 && continue
            site1 = l[i]
            for j in 1:i - 1
                occj = sum(@view occ[(j - 1) * N + 1: j * N])
                occj == 0 && continue
                site2 = translate_to_nearest(l, site1, l[j])
                int_energy += f(site1, site2)::Number * occi * occj
            end
            int_energy += f(site1, site1)::Number * occi * (occi - 1) / 2
        end
        push!(diags, int_energy)
    end
    diagonaloperator(basis(sys), diags)
end

"""
    interaction(f, [T, ]sys, K[; affect_internal=true])

Create an `2K`-site interaction operator for a given `NParticles` system. The function `f`
takes two `K`-tuples of integer numbers, which are site indices for creation and annihilation
operators, and returns the interaction energy.

If `affect_internal` is `true` (default), the interaction operator will act on the internal
degrees of freedom as well, and `f` will take four `K`-tuples - lattice and internal indices for
creation and annihilation operators. If the system has no internal degrees of freedom,
`affect_internal` will automatically be set to `false` and `f` will take two `K`-tuples.
"""
function interaction(f::Function, T::Type{<:Number}, sys::NParticles, ::Val{K}; affect_internal=true) where K
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
    buffer = QuantumOpticsBase.allocate_buffer(bas.occupations)
    for a_inds_c in c_indices
        a_inds = Tuple(a_inds_c)
        issorted(a_inds) || continue
        for at_inds_c in c_indices
            at_inds = Tuple(at_inds_c)
            issorted(at_inds) || continue
            for (m, occ) in enumerate(bas.occupations)
                C = QuantumOpticsBase.state_transition!(buffer, occ, at_inds, a_inds)
                C === nothing && continue
                value = if hasinternal(sys) && affect_internal
                    f(l_inds(a_inds), i_inds(a_inds), l_inds(at_inds), i_inds(at_inds))::Number
                elseif hasinternal(sys)
                    (i_inds(a_inds) == i_inds(at_inds)) * f(l_inds(a_inds), l_inds(at_inds))::Number
                else
                    f(l_inds(a_inds), l_inds(at_inds))::Number
                end
                value == 0 && continue
                n = QuantumOpticsBase.state_index(bas.occupations, buffer)
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
