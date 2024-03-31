using Logging
using ProgressMeter
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs
import KrylovKit: eigsolve

abstract type AbstractEigensystem{BT} end
"""
    Eigensystem{LT, MT} where {LT<:AbstractLattice, MT<:AbstractMatrix}

Eigenvalues and eigenvectors for some operator.
"""
struct Eigensystem{BT<:Basis,MT<:AbstractMatrix} <: AbstractEigensystem{BT}
    basis::BT
    states::MT
    values::Vector{Float64}
    function Eigensystem(basis::BT, states::MT, values::AbstractVector) where {BT,MT}
        @check_size states (length(basis), length(values))
        if issorted(values, by=real)
            new{BT,MT}(basis, states, values)
        else
            sp = sortperm(values, by=real)
            new{BT,MT}(basis, states[:, sp], values[sp])
        end
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
HamiltonianEigensystem(sys::System, eig::AbstractEigensystem) =
    HamiltonianEigensystem(sys, eig.basis, eig.states, eig.values)
function HamiltonianEigensystem(sys::System)
    bas = basis(sys)
    states = zeros(ComplexF64, length(bas), 0)
    values = zeros(ComplexF64, 0)
    return HamiltonianEigensystem(sys, bas, states, values)
end
Eigensystem(ham_eig::HamiltonianEigensystem) =
    Eigensystem(ham_eig.basis, ham_eig.states, ham_eig.values)

function diagonalize_routine(op::DataOperator, ::Val{:lapack}; warning=true, v0=nothing, n=nothing)
    v = size(op.data, 1)
    v > 10000 && warning &&
        @warn """$v×$v dense operator is too large for exact diagonalization; consider making is sparse with `sparse(op)`.
    Set `warning=false` to disable this warning."""
    vals, vecs = eigen(op.data)
    Eigensystem(basis(op), vecs, vals)
end
function diagonalize_routine(op::DataOperator, ::Val{:krylovkit};
        n=10, v0=rand(ComplexF64, size(op.data, 1)), warning=true, kw...)
    v = size(op.data, 1)
    v < 1000 && warning &&
        @warn """$v×$v sparse operator can be diagonalized exactly; consider making is dense with `dense(op)`.
    Set `warning=false` to disable this warning."""
    vals, vecs = eigsolve(op.data, v0, n, :SR, kw...)
    Eigensystem(basis(op), hcat(vecs...), vals)
end
diagonalize_routine(op::DataOperator, ::Val{:auto}; kw...) = diagonalize_routine(find_routine(op)...; kw...)
diagonalize_routine(::DataOperator, ::Val{Routine}) where Routine =
    error("Unsupported diagonalization routine $Routine")

"""
    diagonalize(op::DataOperator[, routine; params...])

Finds eigenvalues and eigenvectors for a `Operator` and stores them in an `Eigensystem`.

Two routines are available:
- `:lapack` uses the `eigen` function from the standard `LinearAlgebra` package.
- `:krylovkit` uses the Lanczos algorithm from the `KrylovKit` package.
    Accepts following parameters:
    - `v0` is the starting vector. Default is `rand(ComplexF64, size(op.data, 1))`.
    - `n` is the target number of eigenvectors. Default is 10.
    All other keyword arguments are passed to the `KrylovKit.eigsolve` function. See its documentation for details.
- `:auto` automatically selects the routine based on the size of the operator.

The default routine is `:lapack` for dense operators. If the operator matrix is less than
5000×5000, it is automatically converted to a dense operator. In other cases `:krylovkit`
is used.

## Example
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(4, 4);

julia> H = tightbinding_hamiltonian(l)
Hamiltonian(dim=16x16)
System: One particle on 16-site 2-dim Bravais lattice in 2D space
16×16 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 48 stored entries:
⎡⠪⡢⠑⢄⠀⠀⠀⠀⎤
⎢⠑⢄⠪⡢⠑⢄⠀⠀⎥
⎢⠀⠀⠑⢄⠪⡢⠑⢄⎥
⎣⠀⠀⠀⠀⠑⢄⠪⡢⎦

julia> eig = diagonalize(H)
Diagonalized Hamiltonian (16 eigenvectors)
Eigenvalues in range -3.23607 .. 3.23607
System: One particle on 16-site 2-dim Bravais lattice in 2D space
```
"""
function diagonalize(ham::Hamiltonian, routine::Val; kw...)
    eig = diagonalize_routine(Operator(ham), routine; kw...)
    HamiltonianEigensystem(ham.sys, eig.basis, eig.states, eig.values)
end
diagonalize(op::DataOperator, routine::Val; kw...) =
    diagonalize_routine(op, routine; kw...)
diagonalize(op::DataOperator, routine::Symbol; kw...) =
    diagonalize(op, Val(routine); kw...)
diagonalize(op::DataOperator; routine=:auto, kw...) =
    diagonalize(op, routine; kw...)
find_routine(op::DenseOpType) = op, Val(:lapack)
find_routine(op::DataOperator) =
    size(op.data, 1) > 5000 ? (op, Val(:krylovkit)) : (dense(op), Val(:lapack))

Base.length(eig::AbstractEigensystem) = length(eig.values)
Base.getindex(eig::AbstractEigensystem, i::Int) = Ket(eig.basis, eig.states[:, i])
Base.getindex(eig::AbstractEigensystem; value::Number) = Ket(eig.basis, eig.states[:,
        argmin(i -> abs(value - eig.values[i]), eachindex(eig.values))])
Base.getindex(eig::AbstractEigensystem, mask::AbstractVector) =
    Eigensystem(eig.basis, eig.states[:, mask], eig.values[mask])
Base.getindex(eig::HamiltonianEigensystem, mask::AbstractVector) =
    HamiltonianEigensystem(eig.sys, eig.basis, eig.states[:, mask], eig.values[mask])
Base.extrema(eig::AbstractEigensystem) = eig.values[begin], eig.values[end]
sample(eig::AbstractEigensystem) = sample(eig.basis)
sample(eig::HamiltonianEigensystem) = sample(eig.sys)
QuantumOpticsBase.basis(eig::AbstractEigensystem) = eig.basis
function Base.show(io::IO, mime::MIME"text/plain", eig::AbstractEigensystem)
    if eig isa Eigensystem
        print(io, "Eigensystem")
    elseif eig isa HamiltonianEigensystem
        print(io, "Diagonalized Hamiltonian")
    end
    println(io, " (", fmtnum(eig, "eigenvector"), ")")
    requires_compact(io) && return
    print(io, "Eigenvalues in range ", @sprintf("%.5f", minimum(eig.values)),
    " .. ", @sprintf("%.5f", maximum(eig.values)))
    if eig isa HamiltonianEigensystem
        print(io, "\nSystem: ")
        show(io, mime, eig.sys)
    end
end

function orthogonalize(eig::AbstractEigensystem, vectors, tol)
    orth = vectors * vectors' * eig.states
    axpby!(1, eig.states, -1, orth)
    norms = [norm(@view orth[:, i]) for i in 1:size(orth, 2)]
    orth_filter = norms .> tol
    if all(orth_filter)
        orth ./= norms'
        return Eigensystem(eig.basis, orth, eig.values)
    else
        normalized_eig = (orth ./ norms')[:, orth_filter]
        return Eigensystem(eig.basis, normalized_eig, eig.values[orth_filter])
    end
end

function _union(eig1::AbstractEigensystem, eigs::AbstractEigensystem...)
    # This function is used to merge eigensystems. No checks are performed.
    newvals = vcat(eig1.values, (eig.values for eig in eigs)...)
    sp = sortperm(newvals)
    newvecs = hcat(eig1.states, (eig.states for eig in eigs)...)[:, sp]
    return Eigensystem(basis(eig1), newvecs, newvals)
end

function Base.union(eig1::AbstractEigensystem, eig2::AbstractEigensystem; tol=1e-10)
    QuantumOpticsBase.check_samebases(basis(eig1), basis(eig2))
    neig2 = orthogonalize(eig2, eig1.states, tol)
    return _union(eig1, neig2)
end
Base.union(heig1::HamiltonianEigensystem, heig2::HamiltonianEigensystem; kw...) =
    HamiltonianEigensystem(heig1.sys, union(Eigensystem(heig1), Eigensystem(heig2); kw...))

"""
    findgroundstate(eig::HamiltonianEigensystem)
    findgroundstate(ham::Hamiltonian)

Finds the ground state of a Hamiltonian. Returns the energy and the state.

## Example
```julia
eig = diagonalize(ham)
E, ψ = findgroundstate(eig)
```
"""
findgroundstate(eig::HamiltonianEigensystem) = eig.values[1], eig[1]
findgroundstate(ham::Hamiltonian) = findgroundstate(diagonalize(ham))

"""
    groundstate(eig::HamiltonianEigensystem)
    groundstate(ham::Hamiltonian)

Finds the ground state of a Hamiltonian. Returns the state.

## Example
```julia
eig = diagonalize(ham)
ψ = groundstate(eig)
```
"""
groundstate(any) = findgroundstate(any)[2]

"""
    projector(eig::Eigensystem)

returns an `Operator` that projects onto the eigenvectors of the spectrum, defined by the formula below.

``\\hat{\\mathcal{P}} = \\sum_i |\\psi_i⟩⟨\\psi_i|``
"""
projector(eig::AbstractEigensystem) = Operator(eig.basis, eig.states * eig.states')

"""
    projector(f, eig::Eigensystem)

Returns an `Operator` representing a function applied to the diagonalized operator defined by the formula below:

``\\hat{\\mathcal{P}} = \\sum_i f(A_i) |\\psi_i⟩⟨\\psi_i|``
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

const TEMP_TOL = 1e-10
const DISABLE_INFO = "; set `info=false` to disable this message"
function groundstate_densitymatrix(eig::AbstractEigensystem; info=true)
    info && @info "Creating density matrix: ground state" * DISABLE_INFO
    ψ = groundstate(eig)
    return (ψ ⊗ ψ')
end
function gibbs_densitymatrix(eig::AbstractEigensystem; T::Real=0, info=true)
    info && @info "Creating density matrix: Gibbs distribution, T = $T" * DISABLE_INFO
    if abs(T) < TEMP_TOL
        return groundstate_densitymatrix(eig, info=false)
    else
        Z = sum(E -> exp(-E / T), eig.values)
        return projector(E -> exp(-E / T) / Z, eig)
    end
end
function ensemble_densitymatrix(eig::AbstractEigensystem;
        μ::Real=0, mu::Real=μ, T::Real=0, statistics::ParticleStatistics=FermiDirac, info=true)
    info && @info "Creating density matrix: $statistics distribution, T = $T, μ = $mu" * DISABLE_INFO
    return projector(densfun(T, mu, statistics), eig)
end
function fermisphere_densitymatrix(eig::AbstractEigensystem; N::Int, info=true)
    info && @info "Creating density matrix: Fermi sphere, N = $N" * DISABLE_INFO
    length(eig) < N &&
        throw(ArgumentError("cannot build Fermi sphere with $N particles: only $(length(Es)) bands present"))
    length(eig) == N && return projector(eig)
    eig.values[N] ≈ eig.values[N+1] && @warn "degenerate levels on the Fermi sphere"
    return projector(eig[1:N])
end
function fixn_densitymatrix(eig::AbstractEigensystem;
        T::Real=0, N::Int, statistics::ParticleStatistics=FermiDirac, maxiter=1000, atol=√eps(), info=true)
    if abs(T) < TEMP_TOL
        if statistics == BoseEinstein
            # condensate
            return groundstate_densitymatrix(eig, info=info) * N
        else
            # Fermi sphere
            return fermisphere_densitymatrix(eig; N=N, info=info)
        end
    end
    info && @info "Creating density matrix: $statistics distribution, N = $N (μ found automatically), T = $T"  * DISABLE_INFO
    Es = eig.values
    newmu = 0.
    for _ in 1:maxiter
        N_ev = sum(densfun(T, newmu, statistics), Es)
        dN_ev = sum(ddensfun(T, newmu, statistics), Es)
        dMu = (N_ev - N) / dN_ev
        abs(dMu) < atol && return ensemble_densitymatrix(eig; T=T, statistics=statistics, mu=newmu, info=false)
        newmu += dMu
    end
    throw(ArgumentError("did not converge"))
end

"""
    densitymatrix(eig::Eigensystem[; T=0, μ, N, statistics, info=true])

Creates an `Operator` representing a equilibrium density matrix, given the eigensystem `eig`
of the Hamiltonian.

The resulting distribution will be Fermi-Dirac or Bose-Einstein if the `statistics` is
specified, otherwise the Gibbs distribution will be used.

## Keyword arguments
- `T` is the temperature of the system. Default is zero.
- `μ` is the chemical potential. Use `mu` as a synonym if Unicode input is not available.
- `N` is the number of particles. If specified, the chemical potential is found automatically.
- `statistics` defines the particle statistics, either `FermiDirac` or `BoseEinstein`.
- `info` is a boolean flag to enable/disable logging. Default is `true`.

Note that if `eig` is a diagonalized `Hamiltonian`, the `μ`, `N` and `statistics` parameters are inserted automatically.
"""
densitymatrix(ham_eig::HamiltonianEigensystem{<:FixedMu}; kw...) =
    ensemble_densitymatrix(Eigensystem(ham_eig);
        T=ham_eig.sys.T, statistics=ham_eig.sys.statistics, mu=ham_eig.sys.chempotential, kw...)
function densitymatrix(ham_eig::HamiltonianEigensystem{<:FixedN}; kw...)
    N = get(kw, :N, ham_eig.sys.nparticles)
    T = get(kw, :T, ham_eig.sys.T)
    statistics = get(kw, :statistics, ham_eig.sys.statistics)
    return fixn_densitymatrix(ham_eig; N=N, T=T, statistics=statistics, kw...)
end
function densitymatrix(ham_eig::HamiltonianEigensystem{<:OneParticleSystem}; kw...)
    if :N in keys(kw)
        (:μ in keys(kw) || :mu in keys(kw)) && throw(ArgumentError("cannot specify both N and μ"))
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
    densitymatrix(HamiltonianEigensystem(OneParticleSystem(sample(eig)), eig); kw...)
densitymatrix(eig::Eigensystem; T::Real = 0) = gibbs_densitymatrix(eig; T = T)
densitymatrix(ham::DataOperator; kw...) = densitymatrix(diagonalize(ham); kw...)
