import QuantumOpticsBase: basis, check_samebases, manybodyoperator

struct Sample{LT, BasisT}
    latt::LT
    internal::BasisT
end
function Sample(l::AbstractLattice, internal::IT=nothing) where {IT<:Nullable{Basis}}
    Sample{typeof(l), IT}(l, internal)
end
const SampleWithoutInternal{LT} = Sample{LT, Nothing}
const SampleWithInternal{LT} = Sample{LT, <:Basis}

Base.length(sample::Sample) = length(sample.latt) * length(sample.internal)
Base.length(sample::SampleWithoutInternal) = length(sample.latt)
QuantumOpticsBase.basis(sample::SampleWithInternal) = sample.internal ⊗ LatticeBasis(sample.latt)
QuantumOpticsBase.basis(sample::SampleWithoutInternal) = LatticeBasis(sample.latt)

function Base.show(io::IO, mime::MIME"text/plain", sample::Sample)
    hasinternal(sample) && print(io, "(")
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
sample(lb::LatticeBasis) = Sample(lb.lat)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].lat, b.bases[1])
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
        print(io, fmtnum(sys.nparticles, "non-interacting particles"), " on ")
    elseif sys isa FixedMu
        print(io, "Non-interactng particles with fixed μ=",
            trunc(sys.chempotential, digits=2), " on ")
    end
    show(io, mime, sys.sample)
end

abstract type ManyBodySystem{SampleT} <: System{SampleT} end
function manybodyoperator(sys::ManyBodySystem, op::AbstractOperator)
    bas = basis(op)
    sys_bas = onebodybasis(sys)
    if bas == onebodybasis(sys)
        return manybodyoperator(ManyBodyBasis(sys_bas, occupations(sys)), op)
    elseif hasinternal(sys)
        if bas == internal_basis(sys)
            latt_op = one(LatticeBasis(lattice(sys)))
            return manybodyoperator(ManyBodyBasis(sys_bas, occupations(sys)), op ⊗ latt_op)
        elseif bas == LatticeBasis(lattice(sys))
            int_op = internal_one(sys)
            return manybodyoperator(ManyBodyBasis(sys_bas, occupations(sys)), int_op ⊗ op)
        end
    end
    throw(ArgumentError("cannot match basis $bas to lattice or on-site phase space"))
end
function manybodyoperator(sys::ManyBodySystem{<:SampleWithInternal}, op::AbstractMatrix)
    N = internal_length(sys)
    @check_size op (N, N)
    return manybodyoperator(sys, Operator(internal_basis(sys), op))
end

struct NParticles{SampleT} <: ManyBodySystem{SampleT}
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
    noun = sys.statistics == FermiDirac ? "fermion" : "boson"
    print(io, "NParticles(", fmtnum(sys.nparticles, noun), ") on ")
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
        throw(ArgumentError("Both chemical potential `μ` and the particle number `N` are defined"))
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

macro accepts_system(fname, default_basis=nothing)
    esc(quote
        $fname(sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(l::AbstractLattice, bas::Basis, args...; kw...) =
            $fname(Sample(l, bas), args...; kw...)
        $fname(l::AbstractLattice, args...; kw...) =
            $fname(Sample(l, $default_basis), args...; kw...)
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
        $fname(type::Type, l::AbstractLattice, bas::Basis, args...; kw...) =
            $fname(type, Sample(l, bas), args...; kw...)
        $fname(type::Type, l::AbstractLattice, args...; kw...) =
            $fname(type, Sample(l, $default_basis), args...; kw...)
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
Base.:(+)(ham::Hamiltonian{Sys, B}, op::Operator{B, B}) where {Sys, B} = Operator(ham) + op
Base.:(-)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where {Sys, B} = op - Operator(ham)
Base.:(-)(ham::Hamiltonian{Sys, B}, op::Operator{B, B}) where {Sys, B} = Operator(ham) - op

function Base.:(+)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    @assert ham.sys == ham2.sys
    return Hamiltonian(ham.sys, Operator(ham) + Operator(ham2))
end
function Base.:(-)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    @assert ham.sys == ham2.sys
    return Hamiltonian(ham.sys, Operator(ham) - Operator(ham2))
end

QuantumOpticsBase.dense(ham::Hamiltonian) = Hamiltonian(ham.sys, dense(Operator(ham)))
QuantumOpticsBase.sparse(ham::Hamiltonian) = Hamiltonian(ham.sys, sparse(Operator(ham)))

Base.summary(io::IO, ham::Hamiltonian) =
    print(io, "Hamiltonian(dim=", length(ham.basis_l), "x", length(ham.basis_r), ")")
function Base.show(io::IO, mime::MIME"text/plain", ham::Hamiltonian)
    summary(io, ham)
    print(io, "\nSystem: ")
    show(io, mime, ham.sys)
    print(io, "\n")
    show(io, mime, ham.data)
end
