import QuantumOpticsBase: Basis

function interaction(f::Function, T::Type{<:Number}, sys::NParticles)
    l = lattice(sys)
    occups = occupations(sys)
    diags = T[]
    N = length(sys.internal)
    for occ in occups
        int_energy = 0.
        for i in eachindex(l)
            occi = sum(@view occ[(i - 1) * N + 1: i * N])
            occi == 0 && continue
            site1 = l[i]
            for j in 1:i-1
                occj = sum(@view occ[(j - 1) * N + 1: j * N])
                occj == 0 && continue
                site2 = l[j]
                int_energy += f(site1, site2) * occi * occj
            end
            int_energy += f(site1, site1) * occi * (occi - 1) / 2
        end
        push!(diags, int_energy)
    end
    diagonaloperator(basis(sys), diags)
end
interaction(f::Function, sys::NParticles) = interaction(f, ComplexF64, sys)
