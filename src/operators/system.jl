using QuantumOpticsBase: basis, check_samebases

struct Sample{LT, BasisT}
    latt::LT
    internal::BasisT
end
function Sample(l::AbstractLattice, internal::IT=nothing;
        boundaries=nothing) where {IT<:Nullable{Basis}}
    new_l = add_boundaries(l, boundaries)
    Sample{typeof(new_l), IT}(new_l, internal)
end
const SampleWithoutInternal{LT} = Sample{LT, Nothing}
const SampleWithInternal{LT} = Sample{LT, <:Basis}

Base.length(sample::Sample) = length(sample.latt) * length(sample.internal)
Base.length(sample::SampleWithoutInternal) = length(sample.latt)
default_bonds(sample::Sample, arg=Val(1)) = default_bonds(lattice(sample), arg)
QuantumOpticsBase.basis(sample::SampleWithInternal) = sample.internal ⊗ LatticeBasis(sample.latt)
QuantumOpticsBase.basis(sample::SampleWithoutInternal) = LatticeBasis(sample.latt)

function Base.show(io::IO, mime::MIME"text/plain", sample::Sample)
    hasinternal(sample) && print("(")
    summary(io, lattice(sample))
    if hasinternal(sample)
        print(io, ") ⊗ ")
        show(io, mime, internal_basis(sample))
    end
end

"""
Returns the `Sample` of the object.

Define this function for your type to implement `Sample` API.

!!! info
    This function can be considered stable internal API. Feel free to use it in your packages.
"""
sample(lb::LatticeBasis) = Sample(lb.sites)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].sites, b.bases[1])
sample(mb::ManyBodyBasis) = sample(mb.onebodybasis)
sample(state::StateType) = sample(basis(state))
sample(op::AbstractLatticeOperator) = sample(basis(op))
lattice(sample::Sample) = sample.latt
lattice(any) = lattice(sample(any))
internal_basis(sample::SampleWithInternal) = sample.internal
internal_basis(::SampleWithoutInternal) = throw(ArgumentError("Sample has no internal basis"))
internal_basis(any) = internal_basis(sample(any))
internal_length(sample::SampleWithInternal) = length(sample.internal)
internal_length(sample::SampleWithoutInternal) = 1
internal_length(any) = internal_length(sample(any))
internal_one(sample::SampleWithInternal) = one(sample.internal)
internal_one(sample::SampleWithoutInternal) = 1
internal_one(any) = internal_one(sample(any))
hasinternal(::SampleWithInternal) = true
hasinternal(::SampleWithoutInternal) = false
hasinternal(any) = hasinternal(sample(any))

abstract type System{SampleT} end
abstract type OneParticleBasisSystem{SampleT} <: System{SampleT} end
sample(sys::System) = sys.sample
default_bonds(sys::System, arg=Val(1)) = default_bonds(lattice(sys.sample), arg)
struct OneParticleSystem{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    T::Float64
    OneParticleSystem(sample::SampleT, T::Real=0) where SampleT = new{SampleT}(sample, Float64(T))
end
OneParticleSystem(l::AbstractLattice, b::Nullable{Basis}=nothing; T=0) = OneParticleSystem(Sample(l, b), T)

QuantumOpticsBase.tensor(l::AbstractLattice, b::Basis) = OneParticleSystem(l, b)
QuantumOpticsBase.tensor(b::Basis, l::AbstractLattice) = OneParticleSystem(l, b)

@enum ParticleStatistics begin
    FermiDirac = 1
    BoseEinstein = -1
end

struct FixedMu{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    chempotential::Float64
    statistics::ParticleStatistics
    T::Float64
end
function FixedMu(sample::SampleT, μ = nothing, mu = μ; statistics=FermiDirac, T = 0) where SampleT
    FixedMu{SampleT}(sample, mu, statistics, T)
end

struct FixedN{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    nparticles::Int
    statistics::ParticleStatistics
    T::Float64
end
function FixedN(sample::SampleT, N; statistics=FermiDirac, T = 0) where SampleT
    FixedN{SampleT}(sample, N, statistics, T)
end

function Base.show(io::IO, mime::MIME"text/plain", sys::OneParticleBasisSystem)
    if sys isa OneParticleSystem
        print(io, "One particle on ")
    elseif sys isa FixedN
        print(io, "$(sys.nparticles) non-interacting particles on ")
    elseif sys isa FixedMu
        print(io, "Non-interactng particles with fixed μ=",
            trunc(sys.chempotential, digits=2), " on ")
    end
    show(io, mime, sys.sample)
end

struct NParticles{SampleT} <: System{SampleT}
    sample::SampleT
    nparticles::Int
    statistics::ParticleStatistics
    T::Float64
end
function NParticles(sample::SampleT, nparticles; statistics = FermiDirac, T = 0) where SampleT<:Sample
    NParticles{SampleT}(sample, nparticles, statistics, T)
end
NParticles(onep::OneParticleSystem, n; kw...) = NParticles(onep.sample, n;
    T = onep.T, kw...)
NParticles(l::AbstractLattice, n; kw...) = NParticles(Sample(l, nothing), n; kw...)
function Base.show(io::IO, mime::MIME"text/plain", sys::NParticles)
    print(io, sys.nparticles, " interacting ",
        sys.statistics == FermiDirac ? "fermions" : "bosons", " on ")
    show(io, mime, sys.sample)
end

function System(sample::Sample; μ = nothing, mu = μ, N = nothing, T = 0, statistics=FermiDirac)
    if mu !== nothing && N === nothing
        return FixedMu(sample, mu, statistics=statistics, T=T)
    elseif N !== nothing && mu === nothing
        return FixedN(sample, N, statistics=statistics, T=T)
    elseif N === mu === nothing
        return OneParticleSystem(sample, T)
    else
        error("Please define the chemical potential `μ` or the particle number `N` (but not both)")
    end
end
System(onep::OneParticleSystem; kw...) = System(onep.sample; T=onep.T, kw...)
System(args...; μ = nothing, mu = μ, N = nothing, statistics=FermiDirac, T=0, kw...) =
    System(Sample(args...; kw...), mu=mu, N=N, T=T, statistics=statistics)

function occupations(np::NParticles)
    if np.statistics == FermiDirac
        fermionstates(length(np.sample), np.nparticles)
    elseif np.statistics == BoseEinstein
        bosonstates(length(np.sample), np.nparticles)
    end
end

onebodybasis(sys::System) = basis(sys.sample)
QuantumOpticsBase.basis(sys::OneParticleBasisSystem) = onebodybasis(sys)
QuantumOpticsBase.basis(sys::NParticles) = ManyBodyBasis(basis(sys.sample), occupations(sys))

Base.zero(sys::System) = zero(basis(sys))

shift_site(sys::System, site) = shift_site(lattice(sys), site)

macro accepts_system(fname, default_basis=nothing)
    esc(quote
        $fname(sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(l::AbstractLattice, bas::Basis, args...; boundaries=nothing, kw...) =
            $fname(Sample(l, bas, boundaries=boundaries), args...; kw...)
        $fname(l::AbstractLattice, args...; boundaries=nothing, kw...) =
            $fname(Sample(l, $default_basis, boundaries=boundaries), args...; kw...)
    end)
end

macro accepts_t(fname)
    esc(quote
    $fname(args...; kw...) = $fname(ComplexF64, args...; kw...)
    $fname(::Type, ::Type, args...; kw...) =
        throw(MethodError($fname, args))
    end)
end

macro accepts_system_t(fname, default_basis=nothing)
    esc(quote
        $fname(type::Type, sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(type, System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(type::Type, l::AbstractLattice, bas::Basis, args...; boundaries=nothing, kw...) =
            $fname(type, Sample(l, bas, boundaries=boundaries), args...; kw...)
        $fname(type::Type, l::AbstractLattice, args...; boundaries=nothing, kw...) =
            $fname(type, Sample(l, $default_basis, boundaries=boundaries), args...; kw...)
        $fname(args...; kw...) = $fname(ComplexF64, args...; kw...)
        $fname(::Type, ::Type, args...; kw...) =
            throw(MethodError($fname, args))
    end)
end

struct Hamiltonian{SystemT, BasisT, T} <: DataOperator{BasisT, BasisT}
    sys::SystemT
    basis_l::BasisT
    basis_r::BasisT
    data::T
end
function Hamiltonian(sys::System, op::Operator)
    return Hamiltonian(sys, basis(op), basis(op), op.data)
end
QuantumOpticsBase.Operator(ham::Hamiltonian) = Operator(ham.basis_l, ham.data)
sample(ham::Hamiltonian) = sample(ham.sys)

Base.:(*)(op::Operator{B1, B2}, ham::Hamiltonian{Sys, B2}) where {Sys, B1, B2} = op * Operator(ham)
Base.:(*)(ham::Hamiltonian{Sys, B2}, op::Operator{B1, B2}) where {Sys, B1, B2} = Operator(ham) * op
Base.:(+)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where {Sys, B} = op + Operator(ham)
Base.:(+)(ham::Hamiltonian{Sys, B}, op::DataOperator{B, B}) where {Sys, B} = Operator(ham) + op
Base.:(-)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where {Sys, B} = op - Operator(ham)
Base.:(-)(ham::Hamiltonian{Sys, B}, op::DataOperator{B, B}) where {Sys, B} = Operator(ham) - op

function Base.:(+)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    @assert ham.sys == ham2.sys
    return Hamiltonian(ham.sys, Operator(ham) + Operator(ham2))
end
function Base.:(-)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    @assert ham.sys == ham2.sys
    return Hamiltonian(ham.sys, Operator(ham) - Operator(ham2))
end
