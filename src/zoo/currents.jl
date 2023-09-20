import QuantumOpticsBase: check_samebases

"""
DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents{N, HT, DT} <: AbstractCurrents
    hamiltonian::HT
    density::DT

    """
        DensityCurrents(hamiltonian, density_mat)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `density_mat`.
    """
    function DensityCurrents(ham::HT, dens::DT) where {HT<:AbstractLatticeOperator, DT<:AbstractLatticeOperator}
        check_samebases(ham, dens)
        new{internal_length(ham), HT, DT}(ham, dens)
    end
end

sample(curr::DensityCurrents) = sample(curr.hamiltonian)
function Base.getindex(curr::DensityCurrents{N}, i::Int, j::Int) where N
    outp_cur = 0.
    for i′ in Base.OneTo(N), j′ in Base.OneTo(N)
        i′′ = i′ + (i - 1) * N
        j′′ = j′ + (j - 1) * N
        outp_cur -= 2 * imag(curr.density.data[i′′, j′′] * curr.hamiltonian.data[j′′, i′′])
    end
    return outp_cur
end
