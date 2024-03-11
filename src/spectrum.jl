using Logging
using ProgressMeter
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
Base.getindex(eig::AbstractEigensystem; value::Number) =
    Ket(eig.basis, eig.states[:, argmin(@. abs(value - eig.values))])
Base.getindex(eig::AbstractEigensystem, mask::AbstractVector) =
    Eigensystem(eig.basis, eig.states[:, mask], eig.values[mask])
Base.getindex(eig::HamiltonianEigensystem, mask::AbstractVector) =
    HamiltonianEigensystem(eig.sys, eig.basis, eig.states[:, mask], eig.values[mask])
Base.extrema(eig::AbstractEigensystem) = eig.values[begin], eig.values[end]
sample(eig::AbstractEigensystem) = sample(eig.basis)
sample(eig::HamiltonianEigensystem) = sample(eig.sys)
QuantumOpticsBase.basis(eig::AbstractEigensystem) = eig.basis

function Base.show(io::IO, ::MIME"text/plain", eig::Eigensystem)
    println(io, "Eigensystem (", fmtnum(eig, "eigenvector"), ")")
    requires_compact(io) && return
    print(io, "Eigenvalues in range $(minimum(eig.values)) .. $(maximum(eig.values))")
end
function Base.show(io::IO, mime::MIME"text/plain", eig::HamiltonianEigensystem)
    println(io, "Diagonalized hamiltonian (", fmtnum(eig, "eigenvector"), ")")
    requires_compact(io) && return
    print(io, "Energies in range $(minimum(eig.values)) .. $(maximum(eig.values))\nSystem: ")
    show(io, mime, eig.sys)
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

const DISABLE_INFO = "; set `info=false` to disable this message"
function groundstate_densitymatrix(eig::AbstractEigensystem; info=true)
    info && @info "Creating density matrix: ground state" * DISABLE_INFO
    ψ = groundstate(eig)
    return (ψ ⊗ ψ')
end
function gibbs_densitymatrix(eig::AbstractEigensystem; T::Real=0, info=true)
    info && @info "Creating density matrix: Gibbs distribution, T = $T" * DISABLE_INFO
    if T == 0
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
    if T ≈ 0
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
        T=ham_eig.sys.T, statistics=statistics, mu=ham_eig.sys.chempotential, kw...)
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

struct GreenFunctionPoint{VecC, VecE}
    weights_up::VecC
    energies_up::VecE
    weights_down::VecC
    energies_down::VecE
    function GreenFunctionPoint(weights_up::VecC, energies_up::VecE,
            weights_down::VecC, energies_down::VecE) where {VecC, VecE}
        @check_size weights_up length(energies_up)
        @check_size weights_down length(energies_down)
        new{VecC, VecE}(weights_up, energies_up, weights_down, energies_down)
    end
end
function (gf::GreenFunctionPoint)(ω::Number)
    return -sum(Base.broadcasted((w, e) -> w / (ω + e), gf.weights_down, gf.energies_down),
        init=-sum(Base.broadcasted((w, e) -> w / (ω - e), gf.weights_up, gf.energies_up)))
end

function Base.show(io::IO, ::MIME"text/plain", gf::GreenFunctionPoint)
    print(io, "Green's function point with ", length(gf.energies_up), " create-bands and ",
        length(gf.energies_down), " annihilate-bands")
end

"""
    GreenFunction

A Green's function for a given lattice and Hamiltonian.
"""
struct GreenFunction{ST, VecC, VecE}
    sample::ST
    weights_up::Matrix{VecC}
    energies_up::VecE
    weights_down::Matrix{VecC}
    energies_down::VecE
    function GreenFunction(sample::ST, mvps::Matrix{VecT}, eps::VecE, mvms=nothing, ems=nothing; E₀=0, E0=E₀) where
            {ST<:Sample, VecT<:AbstractVector, VecE<:AbstractVector}
        p = length(sample)
        @check_size mvps (p,p)
        if mvms === nothing && ems === nothing
            mvms = fill(empty(first(mvps)), p, p)
            ems = empty(eps)
        else
            @check_size mvms (p,p)
        end
        new{ST, VecT, VecE}(sample, mvps, eps .- E0, mvms, ems .- E0)
    end
    function GreenFunction(mvps::Matrix{VecT}, eps::VecE, mvms=nothing, ems=nothing) where
        {VecT<:AbstractVector, VecE<:AbstractVector}
        @check_size mvps :square
        p = size(mvps, 1)
        if mvms === nothing && ems === nothing
            mvms = fill(empty(first(mvps)), p, p)
            ems = empty(eps)
        else
            @check_size mvms (p, p)
        end
        new{Nothing, VecT, VecE}(nothing, mvps, eps, mvms, ems)
    end
end
sample(gf::GreenFunction) = gf.sample
sample(::GreenFunction{Nothing}) = throw(ArgumentError("GreenFunction has no lattice defined"))

function Base.getindex(gf::GreenFunction, α::Int, β::Int)
    return GreenFunctionPoint(gf.weights_up[α, β], gf.energies_up,
        gf.weights_down[α, β], gf.energies_down)
end
function Base.getindex(gf::GreenFunction, site1::AbstractSite, site2::AbstractSite)
    l = lattice(gf)
    i1 = site_index(l, site1)
    i1 === nothing && throw(ArgumentError("site1 is not in the lattice"))
    i2 = site_index(l, site2)
    i2 === nothing && throw(ArgumentError("site2 is not in the lattice"))
    if hasinternal(gf)
        N = internal_length(gf)
        is1 = (i1 - 1) * N + 1:i1 * N
        is2 = (i2 - 1) * N + 1:i2 * N
        return GreenFunction(gf.weights_up[is1, is2], gf.energies_up,
            gf.weights_down[is1, is2], gf.energies_down)
    else
        return gf[i1, i2]
    end
end
inflate_inds(is, N) = N == 1 ? is : [(i - 1) * N + j for i in is for j in 1:N]
function Base.getindex(gf::GreenFunction, any)
    l = lattice(gf)
    inds = to_inds(l, any)
    new_sample = sample(gf)[inds]
    N = internal_length(gf)
    inflated_inds = inflate_inds(inds, N)
    return GreenFunction(new_sample, gf.weights_up[inflated_inds, inflated_inds], gf.energies_up,
        gf.weights_down[inflated_inds, inflated_inds], gf.energies_down)
end
function (gf::GreenFunction)(ω::Number)
    le = length(sample(gf))
    return GreenFunctionEval(gf.sample, [gf[α, β](ω) for α in 1:le, β in 1:le])
end
function (gf::GreenFunction{Nothing})(ω::Number)
    le = size(gf.weights_up, 1)
    return [gf[α, β](ω) for α in 1:le, β in 1:le]
end

function Base.show(io::IO, mime::MIME"text/plain", gf::GreenFunction)
    io = IOContext(io, :compact => true)
    print(io, "Green's function for ")
    show(io, mime, lattice(gf))
end

struct GreenFunctionEval{ST, MT}
    sample::ST
    values::MT
    function GreenFunctionEval(sample::ST, values::MT) where {ST<:Sample,MT}
        @check_size values (length(sample), length(sample))
        new{ST, MT}(sample, values)
    end
end
sample(gf::GreenFunctionEval) = gf.sample
Base.getindex(gf::GreenFunctionEval, α::Int, β::Int) = gf.values[α, β]
function Base.getindex(gf::GreenFunctionEval, site1::AbstractSite, site2::AbstractSite)
    l = lattice(gf)
    i1 = site_index(l, site1)
    i1 === nothing && throw(ArgumentError("site1 is not in the lattice"))
    i2 = site_index(l, site2)
    i2 === nothing && throw(ArgumentError("site2 is not in the lattice"))
    N = internal_length(gf)
    if hasinternal(gf)
        return gf.values[(i1-1)*N+1:i1*N, (i2-1)*N+1:i2*N]
    else
        return gf.values[i1, i2]
    end
end
QuantumOpticsBase.Operator(gf::GreenFunctionEval) = Operator(basis(sample(gf)), gf.values)

"""
    diagonalelements(gf::GreenFunctionEval)

Return the diagonal elements of the Green's function as a `LatticeValue`.
"""
diagonalelements(gf::GreenFunctionEval{<:System{<:SampleWithoutInternal}}) =
    LatticeValue(lattice(gf), diag(gf.values))

function Base.show(io::IO, mime::MIME"text/plain", gf::GreenFunctionEval)
    print(io, "Evaluated Green's function for ")
    summary(io2, lattice(gf))
    requires_compact(io) && return
    print(io, "Values in a ")
    show(io, mime, gf.values)
end

"""
    greenfunction(ham_eig::HamiltonianEigensystem)

Creates a Green's function for a given one-body Hamiltonian eigensystem.
"""
function greenfunction(l, hameig::HamiltonianEigensystem{<:OneParticleBasisSystem}; E₀=0, E0=E₀) # Ensemble formula
    inds = to_inds(lattice(hameig), l)
    N = internal_length(hameig)
    inflated_inds = inflate_inds(inds, N)
    basis(hameig) isa OneParticleBasis ||
        throw(ArgumentError("HamiltonianEigensystem must be on a one-particle basis"))
    GreenFunction(sample(hameig.sys)[inds],
        [conj.(@view hameig.states[α, :]) .* (@view hameig.states[β, :]) for α in inflated_inds, β in inflated_inds],
        hameig.values; E0=E0)
end
function greenfunction(hameig::HamiltonianEigensystem{<:OneParticleBasisSystem}; E₀=0, E0=E₀)
    greenfunction(lattice(hameig), hameig; E0=E0)
end

function _to(state, index, bas::ManyBodyBasis; create)
    bf = basis(state)::ManyBodyBasis
    zs = zeros(ComplexF64, length(bas))
    buffer = QuantumOpticsBase.allocate_buffer(bas.occupations)
    for (i, occ) in enumerate(bf.occupations)
        copyto!(buffer, occ)
        buffer[index] += 1
        if create
            C = QuantumOpticsBase.state_transition!(buffer, occ, (), index) # create
        else
            C = QuantumOpticsBase.state_transition!(buffer, occ, index, ()) # destroy
        end
        C === nothing && continue
        j = QuantumOpticsBase.state_index(bas.occupations, buffer)
        if j === nothing
            @warn "Cannot find state in the basis; check particle numbers"
            continue
        end
        zs[j] = C * state.data[i]
    end
    return Ket(bas, zs)
end

"""
    greenfunction(psi0, hamp, hamm[; E₀, tol, kw...])

Calculates the Green's function for a many-body system with a given initial state `psi0`.

## Arguments
- `psi0` is the initial state.
- `hamp` is the Hamiltonian for the subspace with one more particle than in `psi0`.
- `hamm` is the Hamiltonian for the subspace with one less particle than in `psi0`.

## Keyword arguments
- `E₀` is the energy shift for the Green's function. Default is `0`. Use `E0` as a synonym
    if Unicode input is not available.
- `tol` is the tolerance for the new eigenvectors. Default is `1e-5`.
All other keyword arguments are passed to the `diagonalize` function. See its documentation for details.
"""
function greenfunction(l, psi0::Ket, hamp::Hamiltonian, hamm::Hamiltonian;
        routine=:auto, E₀=0, E0=E₀, tol=1e-5, kw...)
    # Checks...
    bas = basis(psi0)
    bas isa ManyBodyBasis || throw(ArgumentError("`psi0` must be on a many-body basis"))
    obb = bas.onebodybasis
    basis(hamp) isa OneParticleBasis && throw(ArgumentError("`hamp` must be on a manybody basis"))
    obb == basis(hamp).onebodybasis || throw(ArgumentError("`hamp` must be on the same one-particle basis as `psi0`"))
    basis(hamm) isa OneParticleBasis && throw(ArgumentError("`hamm` must be on a manybody basis"))
    obb == basis(hamm).onebodybasis || throw(ArgumentError("`hamm` must be on the same one-particle basis as `psi0`"))

    inds = to_inds(lattice(psi0), l)
    N = internal_length(psi0)
    inflated_inds = inflate_inds(inds, N)

    ep = HamiltonianEigensystem(hamp.sys)
    em = HamiltonianEigensystem(hamm.sys)
    psips = typeof(psi0.data)[]
    psims = typeof(psi0.data)[]
    p = Progress(length(inflated_inds), dt=0.25, desc="Computing GreenFunction...",
        barglyphs=BarGlyphs("[=> ]"))
    for i in inflated_inds
        psip = _to(psi0, i, basis(hamp); create=true)
        push!(psips, psip.data)
        if size(ep.states, 2) != length(basis(hamp))
            psip2 = psip.data - ep.states * ep.states' * psip.data
            if norm(psip2) > tol
                newep = diagonalize(hamp, routine; v0=psip2, kw...)
                ep = union(ep, newep)
            end
        end
        psim = _to(psi0, i, basis(hamm); create=false)
        push!(psims, psim.data)
        if size(em.states, 2) != length(basis(hamm))
            psim2 = psim.data - em.states * em.states' * psim.data
            if norm(psim2) > tol
                newem = diagonalize(hamm, v0=psim2; kw...)
                em = union(em, newem)
            end
        end
        next!(p)
    end
    _weights(states, psis) =
        [conj.(states' * psis[α]) .* (states' * psis[β])
        for α in 1:length(psis), β in 1:length(psis)]
    return GreenFunction(sample(hamp.sys)[inds],
        _weights(ep.states,  psips), ep.values,
        _weights(em.states, psims), em.values; E0=E0)
end
function greenfunction(psi0::Ket, hamp::Hamiltonian, hamm::Hamiltonian; kw...)
    greenfunction(lattice(psi0), psi0, hamp, hamm; kw...)
end

"""
    dos(eig[, E; broaden])
    dos(gf[, E; broaden])

Calculates the DOS (density of states) for a given eigensystem at energy `E`.
If `E` is not specified, a function that calculates the DOS at a given energy is returned.

## Arguments
- `eig` is an `Eigensystem` or `HamiltonianEigensystem`.
- `gf` is a `GreenFunction`.
- `E` is the energy at which the DOS is calculated.

## Keyword arguments
- `broaden` is the broadening factor for the energy levels, default is `0.1`.
"""
dos(eig::AbstractEigensystem, E; broaden=0.1) = imag(sum(1 ./ (eig.values .- (E + im * broaden))))
dos(gf::GreenFunction, E; broaden=0.1) = imag(sum(1 ./ (gf.energies_up .- (E + im * broaden))) -
    sum(1 ./ (gf.energies_down .+ (E + im * broaden))))
dos(any; kw...) = E -> dos(any, E; kw...)

"""
    ldos(gf::GreenFunction, E[; broaden])

Calculates the LDOS (local density of states) for a given Green's function at energy `E`.
`broaden` is the broadening factor for the energy levels, default is `0.1`.
"""
function ldos(gf::GreenFunction, E::Real; broaden=0.1)
    le = length(lattice(gf))
    N = internal_length(gf)
    vals = [sum(imag(gf[α, α](E - im * broaden)) for α in (a - 1) * N + 1:a * N) for a in 1:le]
    return LatticeValue(lattice(gf), vals)
end