using QuantumOpticsBase: basis, check_samebases
struct Sample{LT, BasisT}
    latt::LT
    internal::BasisT
end
Sample(l::LT) where {LT<:AbstractLattice} = Sample{LT,Nothing}(l, nothing)
function Sample(l::BravaisLattice, internal::IT=nothing;
        boundaries=l.boundaries) where {IT<:Nullable{Basis}}
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

sample(lb::LatticeBasis) = Sample(lb.sites)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].sites, b.bases[1])
sample(b::Basis) = throw(MethodError(sample, (b,)))
sample(any) = sample(basis(any))
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

struct NParticles{SampleT} <: System{SampleT}
    sample::SampleT
    nparticles::Int
    statistics::ParticleStatistics
    T::Float64
end
function NParticles(sample::SampleT, nparticles; statistics = FermiDirac, T = 0) where SampleT
    NParticles{SampleT}(sample, nparticles, statistics, T)
end
NParticles(onep::OneParticleSystem, nparticles; kw...) = NParticles(onep.sample, nparticles;
    T = onep.T, kw...)

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

function QuantumOpticsBase.manybodyoperator(ps::NParticles, op::AbstractOperator)
    check_samebases(onebodybasis(ps), basis(op))
    return manybodyoperator(ManyBodyBasis(bas, occupations(ps)), oper)
end
function QuantumOpticsBase.manybodyoperator(sys::OneParticleBasisSystem, op::AbstractOperator)
    check_samebases(onebodybasis(sys), basis(op))
    return op
end
QuantumOpticsBase.manybodyoperator(sys::System, mat::AbstractMatrix) =
    manybodyoperator(sys::System, Operator(onebodybasis(sys), mat))

shift_site(sys::System, site) = shift_site(lattice(sys), site)

macro accepts_system(fname, default_basis=nothing)
    esc(quote
        $fname(sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(l::BravaisLattice, bas::Basis, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, bas, boundaries=boundaries), args...; kw...)
        $fname(l::BravaisLattice, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, $default_basis, boundaries=boundaries), args...; kw...)
    end)
end

macro accepts_system_t(fname, default_basis=nothing)
    esc(quote
        $fname(type::Type, sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(type, System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(type::Type, l::BravaisLattice, bas::Basis, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(type, Sample(l, bas, boundaries=boundaries), args...; kw...)
        $fname(type::Type, l::BravaisLattice, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(type, Sample(l, $default_basis, boundaries=boundaries), args...; kw...)
        $fname(args...; kw...) = $fname(ComplexF64, args...; kw...)
        $fname(::Type, ::Type, args...; kw...) =
            throw(MethodError($fname, args))
    end)
end
