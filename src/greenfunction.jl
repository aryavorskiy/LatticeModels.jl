import QuantumOpticsBase: allocate_buffer, state_transition!, state_index

struct GreenFunctionElement{VecC, VecE}
    weights_up_l::VecC
    weights_up_r::VecC
    energies_up::VecE
    weights_down_l::VecC
    weights_down_r::VecC
    energies_down::VecE
    statistics::ParticleStatistics
    function GreenFunctionElement(weights_up_l::VecC, weights_up_r::VecC, energies_up::VecE,
        weights_down_l::VecC, weights_down_r::VecC, energies_down::VecE, statistics) where {VecC, VecE}
        @check_size weights_up_l length(energies_up)
        @check_size weights_up_r length(energies_up)
        @check_size weights_down_l length(energies_down)
        @check_size weights_down_r length(energies_down)
        new{VecC, VecE}(weights_up_l, weights_up_r, energies_up, weights_down_l, weights_down_r, energies_down, statistics)
    end
end
function (gf::GreenFunctionElement)(ω::Number)
    sum_up = sum(Base.broadcasted((wl, wr, e) -> wl' * wr / (ω - e),
        gf.weights_up_l, gf.weights_up_r, gf.energies_up), init=zero(ComplexF64))
    sum_down = sum(Base.broadcasted((wl, wr, e) -> wl' * wr / (ω + e),
        gf.weights_down_l, gf.weights_down_r, gf.energies_down), init=zero(ComplexF64))
    return sum_up + sum_down * Int(gf.statistics)
end

function Base.show(io::IO, ::MIME"text/plain", gf::GreenFunctionElement)
    print(io, "Green's function element (", length(gf.energies_up), " â†, ",
        length(gf.energies_down), " â bands")
    if length(gf.energies_down) > 0
        print(io, ", ", gf.statistics, " statistics)")
    else    # one particle -> no statistics
        print(io, ")")
    end
end

"""
    GreenFunction

A Green's function for a given lattice and Hamiltonian.
"""
struct GreenFunction{ST, VecC, VecE}
    sample::ST
    weights_up_l::Vector{VecC}
    weights_up_r::Vector{VecC}
    energies_up::VecE
    weights_down_l::Vector{VecC}
    weights_down_r::Vector{VecC}
    energies_down::VecE
    statistics::ParticleStatistics
    function GreenFunction(sample::ST, mvps_l::Vector{VecT}, mvps_r::Vector{VecT}, eps::VecE,
        mvms_l=nothing, mvms_r=nothing, ems=nothing; E₀=0, E0=E₀, statistics) where
            {ST<:Sample, VecT<:AbstractVector, VecE<:AbstractVector}
        p = length(sample)
        @check_size mvps_l p
        @check_size mvps_r p
        if mvms_l === nothing && mvms_r === nothing && ems === nothing
            mvms_l = fill(empty(first(mvps_l)), p)
            mvms_r = fill(empty(first(mvps_r)), p)
            ems = empty(eps)
        elseif mvms_l === nothing || mvms_r === nothing || ems === nothing
            throw(ArgumentError("mvms_l, mvms_r, and ems must be either all or none specified"))
        else
            @check_size mvms_l p
            @check_size mvms_r p
        end
        new{ST, VecT, VecE}(sample, mvps_l, mvps_r, eps .- E0, mvms_l, mvms_r, ems .- E0, statistics)
    end
    function GreenFunction(mvps_l::Vector{VecT}, mvps_r::Vector{VecT}, eps::VecE,
        mvms_l::Vector{VecT}, mvms_r::Vector{VecT}, ems::VecE; statistics) where
        {VecT<:AbstractVector, VecE<:AbstractVector}
        p = length(mvps_l)
        @check_size mvps_r p
        @check_size mvms_l p
        @check_size mvms_r p
        new{Nothing, VecT, VecE}(nothing, mvps_l, mvps_r, eps, mvms_l, mvms_r, ems, statistics)
    end
end
sample(gf::GreenFunction) = gf.sample
sample(::GreenFunction{Nothing}) = throw(ArgumentError("GreenFunction has no lattice defined"))

function Base.getindex(gf::GreenFunction, α::Int, β::Int)
    return GreenFunctionElement(gf.weights_up_l[α], gf.weights_up_r[β], gf.energies_up,
        gf.weights_down_l[α], gf.weights_down_r[β], gf.energies_down, gf.statistics)
end
Base.getindex(gf::GreenFunction, ind1::SampleIndex, ind2::SampleIndex) =
    gf[to_index(gf, ind1), to_index(gf, ind2)]
function Base.getindex(gf::GreenFunction{<:SampleWithInternal},
        site1::AbstractSite, site2::AbstractSite)
    l = lattice(gf)
    i1 = site_index(l, site1)
    i1 === nothing && throw(ArgumentError("site1 is not in the lattice"))
    i2 = site_index(l, site2)
    i2 === nothing && throw(ArgumentError("site2 is not in the lattice"))
    N = internal_length(gf)
    is1 = (i1 - 1) * N + 1:i1 * N
    is2 = (i2 - 1) * N + 1:i2 * N
    return GreenFunction(gf.weights_up_l[is1], gf.weights_up_r[is2], gf.energies_up,
        gf.weights_down_l[is1], gf.weights_down_r[is2], gf.energies_down, statistics=gf.statistics)
end
inflate_inds(is, N) = N == 1 ? is : [(i - 1) * N + j for i in is for j in 1:N]
function Base.getindex(gf::GreenFunction, any)
    l = lattice(gf)
    inds = to_inds(l, any)
    new_sample = sample(gf)[inds]
    N = internal_length(gf)
    inflated_inds = inflate_inds(inds, N)
    return GreenFunction(new_sample, gf.weights_up_l[inflated_inds], gf.weights_up_r[inflated_inds],
        gf.energies_up, gf.weights_down_l[inflated_inds], gf.weights_down_r[inflated_inds],
        gf.energies_down, statistics=gf.statistics)
end
function (gf::GreenFunction)(ω::Number)
    le = length(sample(gf))
    return GreenFunctionEval(gf.sample, [gf[α, β](ω) for α in 1:le, β in 1:le])
end
function (gf::GreenFunction{Nothing})(ω::Number)
    le = size(gf.weights_up_l, 1)
    return [gf[α, β](ω) for α in 1:le, β in 1:le]
end

function Base.show(io::IO, mime::MIME"text/plain", gf::GreenFunction)
    io = IOContext(io, :compact => true)
    print(io, "Green's function for ")
    show(io, mime, sample(gf))
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
Base.getindex(gf::GreenFunctionEval, ind1::SampleIndex, ind2::SampleIndex) =
    gf[to_index(gf, ind1), to_index(gf, ind2)]
function Base.getindex(gf::GreenFunctionEval{<:SampleWithInternal},
        site1::AbstractSite, site2::AbstractSite)
    l = lattice(gf)
    i1 = site_index(l, site1)
    i1 === nothing && throw(ArgumentError("site1 is not in the lattice"))
    i2 = site_index(l, site2)
    i2 === nothing && throw(ArgumentError("site2 is not in the lattice"))
    N = internal_length(gf)
    return gf.values[(i1-1)*N+1:i1*N, (i2-1)*N+1:i2*N]
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
    show(io, mime, sample(gf))
    requires_compact(io) && return
    println(io)
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
    vecs = [hameig.states[α, :] for α in inflated_inds]
    GreenFunction(sample(hameig.sys)[inds], vecs, vecs, hameig.values; E0=E0, statistics=hameig.sys.statistics)
end
function greenfunction(hameig::HamiltonianEigensystem{<:OneParticleBasisSystem}; E₀=0, E0=E₀)
    greenfunction(lattice(hameig), hameig; E0=E0)
end

function _to(state, index, bas::ManyBodyBasis; create)
    bf = basis(state)::ManyBodyBasis
    zs = zeros(ComplexF64, length(bas))
    buffer = allocate_buffer(bf)
    for (i, occ) in enumerate(bf.occupations)
        copyto!(buffer, occ)
        if create
            C = state_transition!(buffer, occ, index, ()) # create
        else
            C = state_transition!(buffer, occ, (), index) # destroy
        end
        C === nothing && continue
        j = state_index(bas.occupations, buffer)
        if j === nothing
            @warn "Cannot find state in the basis; check particle numbers"
            continue
        end
        zs[j] = C * state.data[i]
    end
    return Ket(bas, zs)
end

"""
    greenfunction(psi0, hamp, hamm[; E₀, tol, showprogress, kw...])

Calculates the Green's function for a many-body system with a given initial state `psi0`.

## Arguments
- `psi0` is the initial state.
- `hamp` is the Hamiltonian for the subspace with one more particle than in `psi0`.
- `hamm` is the Hamiltonian for the subspace with one less particle than in `psi0`.

## Keyword arguments
- `E₀` is the energy shift for the Green's function. Default is `0`. Use `E0` as a synonym
    if Unicode input is not available.
- `tol` is the tolerance for the new eigenvectors. Default is `1e-5`.
- `showprogress` is a flag to show the progress bar. Default is `true`.
All other keyword arguments are passed to the `diagonalize` function. See its documentation for details.
"""
function greenfunction(l, psi0::Ket, hamp::Hamiltonian, hamm::Hamiltonian;
        routine=:auto, E₀=0, E0=E₀, tol=1e-5, showprogress=true, kw...)
    # Checks...
    bas = basis(psi0)
    bas isa ManyBodyBasis || throw(ArgumentError("`psi0` must be on a many-body basis"))
    obb = bas.onebodybasis
    obb == basis(hamp).onebodybasis ||
        throw(ArgumentError("`hamp` must be on the same one-particle basis as `psi0`"))
    obb == basis(hamm).onebodybasis ||
        throw(ArgumentError("`hamm` must be on the same one-particle basis as `psi0`"))
    sysp = hamp.sys
    sysm = hamm.sys
    if sysp isa NParticles && sysm isa NParticles
        sysp.nparticles == sysm.nparticles + 2 ||
            throw(ArgumentError("""inconsistent particle numbers.
                Expected `Np - Nm == 2`, got $(sysp.nparticles - sysm.nparticles)"""))
    else
        throw(ArgumentError("Hamiltonians must be on a many-body basis"))
    end

    inds = to_inds(lattice(psi0), l)
    N = internal_length(psi0)
    inflated_inds = inflate_inds(inds, N)

    ep = HamiltonianEigensystem(hamp.sys)
    em = HamiltonianEigensystem(hamm.sys)
    psips = typeof(psi0.data)[]
    psims = typeof(psi0.data)[]
    p = Progress(length(inflated_inds), dt=0.25, desc="Computing GreenFunction...",
        barglyphs=BarGlyphs("[=> ]"), enabled=showprogress)
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
    _weights(states, psis) = [states' * psis[α] for α in 1:length(psis)]
    return GreenFunction(sample(hamp.sys)[inds],
        _weights(ep.states, psips), _weights(ep.states, psips), ep.values,
        _weights(em.states, psims), _weights(em.states, psims), em.values;
        E0=E0, statistics=hamp.sys.statistics)
end
function greenfunction(psi0::Ket, hamp::Hamiltonian, hamm::Hamiltonian; kw...)
    greenfunction(lattice(psi0), psi0, hamp, hamm; kw...)
end
function greenfunction(ham::Hamiltonian, hamp::Hamiltonian, hamm::Hamiltonian; kw...)
    E, psi = findgroundstate(ham)
    greenfunction(psi, hamp, hamm; E0 = E, kw...)
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
dos(eig::AbstractEigensystem, E; broaden=0.1) = imag(sum(1 ./ (eig.values .- (E + im * broaden)))) / pi
dos(gf::GreenFunction, E; broaden=0.1) = imag(
    sum(1 ./ (gf.energies_up .- (E + im * broaden))) + Integer(gf.statistics) *
    sum(1 ./ (gf.energies_down .+ (E + im * broaden)))) / pi
dos(any; kw...) = E -> dos(any, E; kw...)

function _ldos(gf::GreenFunction, E::Real, ind::Int; broaden=0.1)
    N = internal_length(gf)
    return sum(imag(gf[α, α](E - im * broaden)) for α in (ind - 1) * N + 1:ind * N) / pi
end
"""
    ldos(gf::GreenFunction, E[, site; broaden])
    ldos(gf::GreenFunction, site[; broaden])

Calculates the LDOS (local density of states) for a given Green's function at energy `E`.
`broaden` is the broadening factor for the energy levels, default is `0.1`.

If `site` is not specified, the LDOS for all sites is calculated and returned as a `LatticeValue`.
Otherwise, the LDOS for the given site is returned as a `Real` value.

If `E` is not specified, a function that calculates the LDOS at site `site` for given energy
is returned.
"""
function ldos(gf::GreenFunction, E::Real, site::AbstractSite; broaden=0.1)
    l = lattice(gf)
    ind = site_index(l, site)
    ind === nothing && throw(ArgumentError("site is not in the lattice"))
    return _ldos(gf, E, ind; broaden=broaden)
end
ldos(gf::GreenFunction, site::AbstractSite; broaden=0.1) = E -> ldos(gf, E, site; broaden=broaden)
function ldos(gf::GreenFunction, E::Real; broaden=0.1)
    vals = [_ldos(gf, E, ind; broaden=broaden) for ind in eachindex(lattice(gf))]
    return LatticeValue(lattice(gf), vals)
end
