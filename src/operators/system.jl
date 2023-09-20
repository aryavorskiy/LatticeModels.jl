using QuantumOpticsBase: basis, check_samebases
struct Sample{LT, BasisT}
    latt::LT
    internal::BasisT
end
Sample(l::LT) where {LT<:AbstractLattice} = Sample{LT,Nothing}(l, nothing)
function Sample(l::BravaisLattice, internal::IT=nothing;
        boundaries=BoundaryConditions()) where {IT<:Nullable{Basis}}
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

sample(lb::LatticeBasis) = Sample(lb.latt)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].latt, b.bases[1])
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
function FixedMu(sample::SampleT, μ; statistics=FermiDirac, T = 0) where SampleT
    FixedMu{SampleT}(sample, μ, statistics, T)
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
function NParticles(sample::SampleT, nparticles; statistics=FermiDirac, T = 0) where SampleT
    NParticles{SampleT}(sample, nparticles, statistics)
end

function System(sample::Sample; μ = nothing, N = nothing, T = 0, statistics=FermiDirac)
    if μ !== nothing && N === nothing
        return FixedMu(sample, μ, statistics=statistics, T=T)
    elseif N !== nothing && μ === nothing
        return FixedN(sample, N, statistics=statistics, T=T)
    else
        error("Please define the chemical potential `μ` or the particle number `N` (but not both)")
    end
end
System(args...; μ = nothing, N = nothing, statistics=FermiDirac, T=0, kw...) =
    System(Sample(args...; kw...), μ=μ, N=N, T=T, statistics=statistics)

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
        $fname(sample::Sample, args...; T = 0, kw...) =
            $fname(OneParticleSystem(sample, T), args...; kw...)
        $fname(l::BravaisLattice, bas::Basis, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, bas, boundaries=boundaries), args...; kw...)
        $fname(l::BravaisLattice, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, $default_basis, boundaries=boundaries), args...; kw...)
    end)
end
