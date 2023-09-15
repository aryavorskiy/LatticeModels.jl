import QuantumOpticsBase: check_samebases

"""
DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents{HT, DT} <: AbstractCurrents
    hamiltonian::HT
    density::DT

    """
        DensityCurrents(hamiltonian, density_mat)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `density_mat`.
    """
    function DensityCurrents(ham::HT, dens::DT) where {HT<:AbstractLatticeOperator, DT<:AbstractLatticeOperator}
        check_samebases(ham, dens)
        new{HT, DT}(ham, dens)
    end
end

function Base.getindex(curr::DensityCurrents, i::Int, j::Int)
    N = internal_length(curr.hamiltonian)
    is = (i - 1) * N + 1: i * N
    js = (j - 1) * N + 1: j * N
    2imag(tr(curr.density.data[is, js] * curr.hamiltonian.data[js, is]))
end
lattice(curr::DensityCurrents) = lattice(curr.hamiltonian)
