import QuantumOpticsBase: Basis

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
                int_energy += f(site1, site2) * occi * occj
            end
            int_energy += f(site1, site1) * occi * (occi - 1) / 2
        end
        push!(diags, int_energy)
    end
    diagonaloperator(basis(sys), diags)
end
interaction(f::Function, sys::NParticles) = interaction(f, ComplexF64, sys)
