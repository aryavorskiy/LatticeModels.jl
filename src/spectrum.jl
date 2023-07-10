using Logging
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs

"""
    Spectrum{LT, MT} where {LT<:Lattice, MT<:AbstractMatrix}

Eigenvalues and eigenvectors for some operator.
"""
struct Spectrum{BT<:Basis,MT<:AbstractMatrix}
    basis::BT
    states::MT
    energies::Vector{Float64}
    function Spectrum(basis::BT, states::MT, energies::AbstractVector) where {BT,MT}
        length(basis) != size(states)[1] && error("inconsistent basis dimensionality")
        length(energies) != size(states)[2] && error("inconsistent energies list length")
        new{BT,MT}(basis, states, energies)
    end
end

const LatticeOperatorMT{MT} = LatticeOperator{BT,<:MT} where {BT}

"""
    spectrum(op::LatticeOperator)

Finds eigenvalues and eigenvectors for a `LatticeOperator` and stores in in a Spectrum.

!!! note
    This method finds eigenvalues and eigenvectors using `LinearAlgebra.eigen`, which can be not defined for some array types.
    Consider redefining it for your array type or constructing the Spectrum object explicitly.
"""
function spectrum(lop::LatticeOperator)
    !all(isfinite.(lop.data)) && error("NaN of Inf in operator matrix")
    vals, vecs = eigen(Hermitian(lop.data))
    Spectrum(basis(lop), vecs, vals)
end

eigvals(sp::Spectrum) = sp.energies
eigvecs(sp::Spectrum) = sp.states
basis(sp::Spectrum) = sp.basis

Base.length(sp::Spectrum) = length(sp.energies)
Base.getindex(sp::Spectrum, i::Int) = LatticeArray(sp.basis, sp.states[:, i])
Base.getindex(sp::Spectrum; E::Number) =
    LatticeArray(sp.basis, sp.states[:, argmin(@. abs(E - sp.energies))])
Base.getindex(sp::Spectrum, mask) =
    Spectrum(sp.basis, sp.states[:, mask], sp.energies[mask])

function Base.show(io::IO, ::MIME"text/plain", sp::Spectrum)
    println(io, "Spectrum with $(length(sp)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(sp.energies)) .. $(maximum(sp.energies))")
end

@doc raw"""
    projector(sp::Spectrum)

$$\hat{\mathcal{P}} = \sum_i |\psi_i⟩⟨\psi_i|$$

Creates a `LatticeOperator` that projects onto the eigenvectors of the spectrum, described by the formula above.
"""
projector(sp::Spectrum) = LatticeArray(sp.basis, sp.states * sp.states')

@doc raw"""
    projector(f, sp::Spectrum)

$$\hat{\mathcal{P}} = \sum_i p_i |\psi_i⟩⟨\psi_i|$$

Creates a `LatticeOperator` that projects onto the eigenvectors of the spectrum, described by the formula above.
The ``p_i`` amplitudes are defined by the `f` function, which takes the eigenvalue ``E_i`` and returns a number (or a boolean).
"""
projector(f::Function, sp::Spectrum) =
    LatticeArray(sp.basis, sp.states * (f.(sp.energies) .* sp.states'))

"""
    filled_projector(sp::Spectrum[, fermi_level=0])

Creates a `LatticeOperator` that projects onto the eigenvectors which have eigenvalues ``E_i`` less than `fermi_level` (0 by default).

Same as `projector(<(fermi_level), sp)`, see
"""
filled_projector(sp::Spectrum, fermi_level=0) = projector(<(fermi_level), sp)

"""
    fermi_dirac(μ, T)

Generates a function that takes the energy and returns the state density acccording to Fermi-Dirac statistics.
"""
fermi_dirac(μ, T) = E -> 1 / (exp((E - μ) / T) + 1)

"""
    bose_einstein(μ, T)

Generates a function that takes the energy and returns the state density acccording to Bose-Einstein statistics.
"""
bose_einstein(μ, T) = E -> 1 / (exp((E - μ) / T) - 1)

@doc raw"""
    dos(sp::Spectrum, δ)

Generates a function to calculate the DOS (Density of States), which is defined as
$\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ and can be understood as a sum
of Lorenz distributions with width equal to $\delta$.
"""
dos(sp::Spectrum, δ::Real) = (E -> imag(sum(1 ./ (eigvals(sp) .- (E + im * δ)))))

@doc raw"""
    ldos(sp::Spectrum, E, δ)

Calculates the LDOS (Local Density of States), which is defined as the imaginary part of partial trace of
$\frac{1}{\hat{H} - E - i\delta}$ operator.
"""
function ldos(sp::Spectrum, E::Real, δ::Real)
    Es = eigvals(sp)
    Vs = eigvecs(sp)
    l = lattice(sp)
    N = dims_internal(sp)
    inves = imag.(1 ./ (Es .- (E + im * δ)))'
    LatticeValue(l, [sum(abs2.(Vs[(i-1)*N+1:i*N, :]) .* inves) for i in 1:length(l)])
end

"""
    ldos(sp::Spectrum, δ)

Generates a function that accepts the energy `E` and returns `ldos(sp, E, δ)`.
Use this if you want to find the LDOS for the same `Spectrum`, but for many different values of `E` -
the produced function is optimized and reduces overall computation time dramatically.
"""
function ldos(sp::Spectrum, δ::Real)
    Es = eigvals(sp)
    l = lattice(sp)
    N = dims_internal(sp)
    density_sums = reshape(
        sum(reshape(abs2.(eigvecs(sp)), (N, :, length(sp))), dims=1), (:, length(sp)))
    E -> LatticeValue(l, vec(sum(density_sums .* imag.(1 ./ (Es .- (E + im * δ)))', dims=2)))
end
