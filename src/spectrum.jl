using Logging
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs
import KrylovKit: eigsolve

abstract type AbstractEigensystem{BT} end
"""
    Eigensystem{LT, MT} where {LT<:Lattice, MT<:AbstractMatrix}

Eigenvalues and eigenvectors for some operator.
"""
struct Eigensystem{BT<:Basis,MT<:AbstractMatrix} <: AbstractEigensystem{BT}
    basis::BT
    states::MT
    values::Vector{Float64}
    function Eigensystem(basis::BT, states::MT, values::AbstractVector) where {BT,MT}
        @check_size states (length(basis), length(values))
        sp = sortperm(real.(values))
        new{BT,MT}(basis, states[:, sp], values[sp])
    end
end
struct HamiltonianEigensystem{ST<:System,BT<:Basis,MT<:AbstractMatrix} <: AbstractEigensystem{BT}
    sys::ST
    basis::BT
    states::MT
    values::Vector{Float64}
    function HamiltonianEigensystem(sys::ST, basis::BT, states::MT, values::AbstractVector) where {ST,BT,MT}
        @check_size states (length(basis), length(values))
        sp = sortperm(real.(values))
        new{ST,BT,MT}(sys, basis, states[:, sp], values[sp])
    end
end
HamiltonianEigensystem(eig::AbstractEigensystem, sys::System) =
    HamiltonianEigensystem(sys, eig.basis, eig.states, eig.values)
Eigensystem(ham_eig::HamiltonianEigensystem) =
    Eigensystem(ham_eig.basis, ham_eig.states, ham_eig.values)
QuantumOpticsBase.basis(eig::AbstractEigensystem) = eig.basis

"""
    diagonalize(op::DataOperator[, routine; params...])

Finds eigenvalues and eigenvectors for a `Operator` and stores them in an `Eigensystem`.

Two routines are available:
- `:lapack` uses the `eigen` function from the standard `LinearAlgebra` package.
- `:krylovkit` uses the Lanczos algorithm from the `KrylovKit` package.
    Accepts following parameters:
    - `v0` is the starting vector. Default is `rand(ComplexF64, size(op.data, 1))`.
    - `n` is the target number of eigenvectors. Default is 10.
    Remaining keyword arguments are passed to the `eigsolve` function.

The default routine is `:lapack` for dense operators. If the operator matrix is less than
3000×3000, it is automatically converted to a dense operator. In other cases `:krylovkit`
is used.
"""
function diagonalize(op::DataOperator, ::Val{:lapack}; warning=true)
    v = size(op.data, 1)
    v > 10000 && warning &&
        @warn """$v×$v dense operator is too large for exact diagonalization; consider making is sparse with `sparse(op)`.
    Set `warning=false` to disable this warning."""
    vals, vecs = eigen(op.data)
    Eigensystem(basis(op), vecs, vals)
end
function diagonalize(op::DataOperator, ::Val{:krylovkit};
        n=10, v0=rand(ComplexF64, size(op.data, 1)), warning=true, kw...)
    v = size(op.data, 1)
    v < 1000 && warning &&
        @warn """$v×$v sparse operator can be diagonalized exactly; consider making is dense with `dense(op)`.
    Set `warning=false` to disable this warning."""
    vals, vecs = eigsolve(op.data, v0, n, :SR, kw...)
    Eigensystem(basis(op), hcat(vecs...), vals)
end
diagonalize(::DataOperator, ::Val{Routine}) where Routine =
    error("Unsupported diagonalization routine $Routine")
diagonalize(op::DataOperator, routine::Symbol; kw...) =
    diagonalize(op, Val(routine); kw...)
function diagonalize(ham::Hamiltonian, routine; kw...)
    eig = diagonalize(Operator(ham), routine; kw...)
    HamiltonianEigensystem(ham.sys, eig.basis, eig.states, eig.values)
end
diagonalize(op::DataOperator; kw...) =
    diagonalize(find_routine(op)...; kw...)
find_routine(op::DenseOpType) = op, Val(:lapack)
find_routine(op::DataOperator) =
    size(op.data, 1) > 3000 ? (op, Val(:krylovkit)) : (dense(op), Val(:lapack))


Base.length(eig::AbstractEigensystem) = length(eig.values)
Base.getindex(eig::AbstractEigensystem, i::Int) = Ket(eig.basis, eig.states[:, i])
Base.getindex(eig::AbstractEigensystem; value::Number) =
    Ket(eig.basis, eig.states[:, argmin(@. abs(value - eig.values))])
Base.getindex(eig::AbstractEigensystem, mask) =
    Eigensystem(eig.basis, eig.states[:, mask], eig.values[mask])
sample(eig::AbstractEigensystem) = sample(eig.basis)

function Base.show(io::IO, ::MIME"text/plain", eig::AbstractEigensystem)
    println(io, "Eigensystem with $(length(eig)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(eig.values)) .. $(maximum(eig.values))")
end

groundstate(eig::HamiltonianEigensystem) = eig[1]

@doc raw"""
    projector(eig::Eigensystem)

returns an `Operator` that projects onto the eigenvectors of the spectrum, defined by the formula below.

$$\hat{\mathcal{P}} = \sum_i |\psi_i⟩⟨\psi_i|$$
"""
projector(eig::AbstractEigensystem) = Operator(eig.basis, eig.states * eig.states')

@doc raw"""
    projector(f, eig::Eigensystem)

Returns an `Operator` representing a function applied to the diagonalized operator defined by the formula below:

$$\hat{\mathcal{P}} = \sum_i f(A_i) |\psi_i⟩⟨\psi_i|$$
"""
function projector(f, eig::AbstractEigensystem)
    Operator(eig.basis, eig.states * (f.(eig.values) .* eig.states'))
end

@inline function densfun(T, mu, statistics::ParticleStatistics)
    return E -> (E - mu == T == 0) ? 1. : 1 / (exp((E - mu) / T) + Int(statistics))
end
@inline function ddensfun(T, mu, statistics::ParticleStatistics)
    return E -> -exp((E - mu) / T) / T / (exp((E - mu) / T) + Int(statistics))^2
end

function groundstate_densitymatrix(eig::AbstractEigensystem)
    ψ = groundstate(eig)
    return (ψ ⊗ ψ')
end
function gibbs_densitymatrix(eig::AbstractEigensystem; T::Real=0)
    if T == 0
        return groundstate_densitymatrix(eig)
    else
        Z = sum(E -> exp(-E / T), eig.values)
        return projector(E -> exp(-E / T) / Z, eig)
    end
end
function ensemble_densitymatrix(eig::AbstractEigensystem;
        μ::Real=0, mu::Real=μ, T::Real=0, statistics::ParticleStatistics=FermiDirac)
    return projector(densfun(T, mu, statistics), eig)
end
function fermisphere_densitymatrix(eig::AbstractEigensystem; N::Int)
    length(Es) < N &&
        throw(ArgumentError("cannot build Fermi sphere with $N particles: only $(length(Es)) bands present"))
    length(Es) == N && return projector(ham_eig)
    Es[N] ≈ Es[N+1] && @warn "degenerate levels on the Fermi sphere"
    return projector(eig[1:N])
end
function fixn_densitymatrix(eig::AbstractEigensystem;
        T::Real=0, N::Int, statistics::ParticleStatistics=FermiDirac, maxiter=1000, atol=√eps())
    Es = eig.values
    T ≈ 0 && @warn "chempotential detection may fail on low temperatures"
    newmu = 0.
    for _ in 1:maxiter
        N_ev = sum(densfun(T, newmu, statistics), Es)
        dN_ev = sum(ddensfun(T, newmu, statistics), Es)
        dMu = (N_ev - N) / dN_ev
        abs(dMu) < atol && return ensemble_densitymatrix(eig; T=T, statistics=statistics, mu=newmu)
        newmu += dMu
    end
    throw(ArgumentError("did not converge"))
end

"""
    densitymatrix(eig::Eigensystem[; T=0, μ=0, statistics])

Creates an `Operator` representing a equilibrium density matrix, given the eigensystem `eig`
of the Hamiltonian.

The resulting distribution will be Fermi-Dirac or Bose-Einstein if the `statistics` is
specified, otherwise the Gibbs distribution will be used.

## Keyword arguments
- `T` is the temperature of the system. Default is zero.
- `μ` is the chemical potential. Use keyword `mu` as a synonym if Unicode input is not available.
- `field` is the magnetic field. Default is `NoField()`.
- `statistics` defines the particle statistics, either `FermiDirac` or `BoseEinstein`.

Note that if `eig` is a diagonalized `Hamiltonian`, the `μ` and `statistics` parameters are inserted automatically.
"""
densitymatrix(ham_eig::HamiltonianEigensystem{<:FixedMu}; kw...) =
    ensemble_densitymatrix(Eigensystem(ham_eig);
        T=ham_eig.sys.T, statistics=statistics, mu=ham_eig.sys.chempotential, kw...)
function densitymatrix(ham_eig::HamiltonianEigensystem{<:FixedN}; kw...)
    N = get(kw, :N, ham_eig.sys.nparticles)
    T = get(kw, :T, ham_eig.sys.T)
    statistics = get(kw, :statistics, ham_eig.sys.statistics)
    if T ≈ 0
        if statistics == BoseEinstein
            # condensate
            return groundstate_densitymatrix(ham_eig) * N
        else
            # Fermi sphere
            return fermisphere_densitymatrix(ham_eig; N=N)
        end
    else
        # Newton's method
        return fixn_densitymatrix(ham_eig; N=N, T=T, statistics=statistics, kw...)
    end
end
function densitymatrix(ham_eig::HamiltonianEigensystem{<:OneParticleSystem}; kw...)
    if :N in keys(kw)
        return fixn_densitymatrix(ham_eig; T = ham_eig.sys.T, kw...)
    elseif :μ in keys(kw) || :mu in keys(kw) || :statistics in keys(kw)
        return ensemble_densitymatrix(ham_eig; T = ham_eig.sys.T, kw...)
    else
        return gibbs_densitymatrix(ham_eig; T = ham_eig.sys.T, kw...)
    end
end
function densitymatrix(ham_eig::HamiltonianEigensystem{<:NParticles}; kw...)
    return gibbs_densitymatrix(ham_eig; T = ham_eig.sys.T, kw...)
end
densitymatrix(eig::Eigensystem{<:AbstractLatticeBasis}; kw...) =
    densitymatrix(HamiltonianEigensystem(eig, OneParticleSystem(sample(eig))); kw...)
densitymatrix(eig::Eigensystem; T::Real = 0) = gibbs_densitymatrix(eig; T = T)
densitymatrix(ham::DataOperator; kw...) = densitymatrix(diagonalize(ham); kw...)

@doc raw"""
    dos(eig::Eigensystem, δ)

Generates a function to calculate the DOS (Density of States), which is defined as
$\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ and can be understood as a sum
of Lorenz distributions with width equal to $\delta$.
"""
dos(eig::AbstractEigensystem, δ::Real) = (E -> imag(sum(1 ./ (eig.values .- (E + im * δ)))))

@doc raw"""
    ldos(eig::Eigensystem, E, δ)

Calculates the LDOS (Local Density of States), which is defined as the imaginary part of partial trace of
$\frac{1}{\hat{H} - E - i\delta}$ operator.
"""
function ldos(eig::AbstractEigensystem{<:AbstractLatticeBasis}, E::Real, δ::Real)
    Es = eig.values
    Vs = eig.states
    l = lattice(eig.basis)
    N = internal_length(eig.basis)
    inves = imag.(1 ./ (Es .- (E + im * δ)))'
    LatticeValue(l, [sum(abs2.(Vs[(i-1)*N+1:i*N, :]) .* inves) for i in eachindex(l)])
end

"""
    ldos(eig::Eigensystem, δ)

Generates a function that accepts the energy `E` and returns `ldos(sp, E, δ)`.
Use this if you want to find the LDOS for the same `Spectrum`, but for many different values of `E` -
the produced function is optimized and reduces overall computation time dramatically.
"""
function ldos(eig::AbstractEigensystem{<:AbstractLatticeBasis}, δ::Real)
    Es = eig.values
    l = lattice(eig.basis)
    N = internal_length(eig.basis)
    density_sums = reshape(
        sum(reshape(abs2.(eig.states), (N, :, length(eig))), dims=1), (:, length(eig)))
    E -> LatticeValue(l, vec(sum(density_sums .* imag.(1 ./ (Es .- (E + im * δ)))', dims=2)))
end
