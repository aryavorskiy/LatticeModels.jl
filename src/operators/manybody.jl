import QuantumOpticsBase: Basis, SparseOperator

function translate_to_nearest(l::BravaisLattice{N}, site1::BravaisSite{N}, site2::BravaisSite{N}) where N
    min_lp = site2.lp
    min_dist = norm(min_lp.unit_cell - site1.lp.unit_cell)
    for cind in cartesian_indices(l)
        tup = Tuple(cind)
        tr_vec = @SVector zeros(Int, N)
        for i in eachindex(tup)
            tr_vec += tup[i] * l.boundaries.bcs[i].R
        end
        new_lp = shift_site(tr_vec, site2.lp)
        if norm(new_lp.unit_cell - site1.lp.unit_cell) < min_dist
            min_dist = norm(new_lp.unit_cell - site1.lp.unit_cell)
            min_lp = new_lp
        end
    end
    return get_site(l, min_lp)
end

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

# This function uses undocumented QuantumOpticsBase API. Be careful!
function interaction(f::Function, T::Type{<:Number}, sys::NParticles, ::Val{K}; affect_internal=true) where K
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
                C = QuantumOpticsBase.state_transition!(buffer, occ, a_inds, at_inds)
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
