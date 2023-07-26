using Logging
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs

"""
    Eigensystem{LT, MT} where {LT<:Lattice, MT<:AbstractMatrix}

Eigenvalues and eigenvectors for some operator.
"""
struct Eigensystem{BT<:Basis,MT<:AbstractMatrix}
    basis::BT
    states::MT
    values::Vector{Float64}
    function Eigensystem(basis::BT, states::MT, values::AbstractVector) where {BT,MT}
        length(basis) != size(states)[1] && error("inconsistent basis dimensionality")
        length(values) != size(states)[2] && error("inconsistent energies list length")
        new{BT,MT}(basis, states, values)
    end
end
QuantumOpticsBase.basis(eig::Eigensystem) = eig.basis

"""
    diagonalize(op::Operator)

Finds eigenvalues and eigenvectors for a `Operator` and stores it in an Eigensystem.
"""
function diagonalize(op::Operator)
    # TODO support sparse diagonalization also
    vals, vecs = eigen(dense(op).data)
    Eigensystem(basis(op), vecs, vals)
end

Base.length(eig::Eigensystem) = length(eig.values)
Base.getindex(eig::Eigensystem, i::Int) = Ket(eig.basis, eig.states[:, i])
Base.getindex(eig::Eigensystem; value::Number) =
    Ket(eig.basis, eig.states[:, argmin(@. abs(value - eig.values))])
Base.getindex(eig::Eigensystem, mask) =
    Eigensystem(eig.basis, eig.states[:, mask], eig.values[mask])

function Base.show(io::IO, ::MIME"text/plain", eig::Eigensystem)
    println(io, "Eigensystem with $(length(eig)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(eig.values)) .. $(maximum(eig.values))")
end

@doc raw"""
    projector(eig::Eigensystem)

returns an `Operator` that projects onto the eigenvectors of the spectrum, defined by the formula below.

$$\hat{\mathcal{P}} = \sum_i |\psi_i⟩⟨\psi_i|$$
"""
projector(eig::Eigensystem) = Operator(eig.basis, eig.states * eig.states')

@doc raw"""
    apply_to_eigenvalues(f, eig::Eigensystem)

Returns an `Operator` representing a function applied to the diagonalized operator defined by the formula below:

$$\hat{\mathcal{P}} = \sum_i f(A_i) |\psi_i⟩⟨\psi_i|$$
"""
function apply_to_eigenvalues(f, eig::Eigensystem)
    Operator(eig.basis, eig.states * (f.(eig.values) .* eig.states'))
end

"""
    densitymatrix(eig::Eigensystem[; μ=0, T=0, statistics=FermiDirac])

Creates an `Operator` representing a equilibrium density matrix, given the eigensystem `eig`
of the Hamiltonian.

## Keyword arguments
The `μ` and `T` keywords set the chemical potential and the temperature respectively.
The `statistics` Keyword sets the probability distribution .`OneParticle`, which is the
default, means Boltzmann distribution.
"""
function densitymatrix(eig::Eigensystem; μ::Real=0, T::Real=0, statistics::ParticleStatistics=OneParticle)
    apply_to_eigenvalues(eig) do E
        E - μ == T == 0 ? 1. : 1 / (exp((E - μ) / T) + Int(statistics))
    end
end

@doc raw"""
    dos(eig::Eigensystem, δ)

Generates a function to calculate the DOS (Density of States), which is defined as
$\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ and can be understood as a sum
of Lorenz distributions with width equal to $\delta$.
"""
dos(eig::Eigensystem, δ::Real) = (E -> imag(sum(1 ./ (eig.values .- (E + im * δ)))))

@doc raw"""
    ldos(eig::Eigensystem, E, δ)

Calculates the LDOS (Local Density of States), which is defined as the imaginary part of partial trace of
$\frac{1}{\hat{H} - E - i\delta}$ operator.
"""
function ldos(eig::Eigensystem{<:AbstractLatticeBasis}, E::Real, δ::Real)
    Es = eig.values
    Vs = eig.states
    l = lattice(eig.basis)
    N = length(internal_basis(eig.basis))
    inves = imag.(1 ./ (Es .- (E + im * δ)))'
    LatticeValue(l, [sum(abs2.(Vs[(i-1)*N+1:i*N, :]) .* inves) for i in 1:length(l)])
end

"""
    ldos(eig::Eigensystem, δ)

Generates a function that accepts the energy `E` and returns `ldos(sp, E, δ)`.
Use this if you want to find the LDOS for the same `Spectrum`, but for many different values of `E` -
the produced function is optimized and reduces overall computation time dramatically.
"""
function ldos(eig::Eigensystem{<:AbstractLatticeBasis}, δ::Real)
    Es = eig.values
    l = lattice(eig.basis)
    N = length(internal_basis(eig.basis))
    density_sums = reshape(
        sum(reshape(abs2.(eig.states), (N, :, length(eig))), dims=1), (:, length(eig)))
    E -> LatticeValue(l, vec(sum(density_sums .* imag.(1 ./ (Es .- (E + im * δ)))', dims=2)))
end
