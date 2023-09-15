import QuantumOpticsBase: Basis, SparseOperator, SparseOpPureType

# Override this function with more performant one
function QuantumOpticsBase.manybodyoperator_1(basis::ManyBodyBasis,
        op::Operator{<:AbstractLatticeBasis, <:AbstractLatticeBasis})
    N = length(basis)
    S = length(basis.onebodybasis)
    builder = SparseMatrixBuilder{ComplexF64}(N, N)
    M = op.data
    @inbounds for m=1:N, n=1:N
        # This code is horrific, but 10x more performant. DO NOT TOUCH!!!
        # Wait until my pull request is accepted, then we can get rid of this
        diffcount = 0
        occ_m = basis.occupations[m]
        occ_n = basis.occupations[n]
        mi = 0
        ni = 0
        for i in 1:S
            occ_m[i] == occ_n[i] && continue
            diffcount += 1
            diffcount ≥ 3 && break
            di = occ_m[i] - occ_n[i]
            if di == 1
                mi = i
            elseif di == -1
                ni = i
            else
                diffcount = 1
                break
            end
        end
        if diffcount == 0
            mat_el = 0.
            for i in 1:S
                occ_m[i] == 0 && continue
                mat_el += occ_m[i] * M[i, i]
            end
            increment!(builder, mat_el, m, n)
        elseif diffcount == 2
            mi != 0 && ni != 0 &&
                increment!(builder,  M[mi, ni], m, n, factor = √(occ_m[mi] * occ_n[ni]))
        end
    end
    return SparseOperator(basis, to_matrix(builder))
end

function interaction(f::Function, T::Type{<:Number}, sample::Sample)
    l = lattice(sample)
    !ismanybody(sample) &&
        throw(ArgumentError("Cannot define interaction on one-particle sample"))
    occups = occupations(sample)
    diags = T[]
    N = length(sample.internal)
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
    diagonaloperator(ManyBodyBasis(basis(sample), occups), diags)
end
interaction(f::Function, sample::Sample) = interaction(f, ComplexF64, sample)
