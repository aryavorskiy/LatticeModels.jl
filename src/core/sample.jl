import QuantumOpticsBase: Basis, AbstractOperator, basis, check_samebases
abstract type Boundary{N} end

shift_site(js::SVector{N, Int}, lp::BravaisPointer{N}) where N =
    BravaisPointer(lp.unit_cell + js, lp.basis_index)
shift_site(js::SVector{N}, site::BravaisSite{N}) where N =
    BravaisSite(shift_site(js, site.lp), site.bravais)

struct TwistedBoundary{N} <: Boundary{N}
    R::SVector{N, Int}
    Θ::Float64
end
PeriodicBoundary(svec) = TwistedBoundary(svec, 0)
function shift_site(bc::TwistedBoundary{N}, i::Int, site) where N
    i == 0 && return 1., site
    return exp(im * i * bc.Θ), shift_site(-bc.R * i, site)
end

struct FunctionBoundary{N, F<:Function} <: Boundary{N}
    condition::F
    R::SVector{N, Int}
end

function shift_site(bc::FunctionBoundary{N}, i::Int, site::BravaisSite{N}) where N
    i == 0 && return 1., site
    for _ in 1:abs(i)
        if i > 0
            site = shift_site(-bc.R, site)
            factor *= bc.condition(site)
        else
            factor /= bc.condition(site)
            site = shift_site(dv * lspan, l, site)
        end
    end
    factor, site
end

struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    function BoundaryConditions(bcs::CondsTuple) where CondsTuple<:NTuple{M, <:Boundary{N}} where {M, N}
        new{CondsTuple}(bcs)
    end
end
BoundaryConditions(args::BoundaryConditions...) = BoundaryConditions(args)

@generated cartesian_indices(depth::Int, ::Val{M}) where M = quote
    CartesianIndex($((:(-depth) for _ in 1:M)...)):CartesianIndex($((:depth for _ in 1:M)...))
end
function route(bcs::BoundaryConditions{<:NTuple{M}}, l::BravaisLattice, lp::BravaisPointer{N}, depth=1) where {M, N}
    for cind in cartesian_indices(depth, Val(M))
        tup = Tuple(cind)
        tr_vec = @SVector zeros(Int, N)
        for i in 1:M
            tr_vec += tup[i] * bcs.bcs[i].R
        end
        new_lp = shift_site(-tr_vec, lp)
        new_lp in l && return tr_vec
    end
    return @SVector zeros(M)
end
function shift_site(bcs::BoundaryConditions, l::BravaisLattice, site)
    factor = 1.
    tr_vec = route(bcs, l, site.lp)
    for i in eachindex(tr_vec)
        new_factor, site = shift_site(bcs.bcs[i], i, site)
        factor *= new_factor
    end
    return factor, site
end

struct MagneticBoundaryConditions end

struct Sample{LT, BasisT, BoundaryT}
    latt::LT
    boundaries::BoundaryT
    internal::BasisT
end
const SampleWithoutInternal{LT, BoundaryT} = Sample{LT, Nothing, BoundaryT}
const SampleWithInternal{LT, BoundaryT} = Sample{LT, <:Basis, BoundaryT}

function Sample(latt::LT, internal::BT=nothing;
        boundaries=BoundaryConditions()) where {LT<:BravaisLattice, BT<:Nullable{Basis}}
    return Sample(latt, boundaries, internal)
end

Base.length(sample::Sample) = length(sample.latt) * length(sample.internal)
Base.length(sample::SampleWithoutInternal) = length(sample.latt)
lattice(sample::Sample) = sample.latt
default_bonds(sample::Sample, arg=Val(1)) = default_bonds(lattice(sample), arg)
internal_one(sample::Sample) = one(sample.internal)
internal_one(sample::SampleWithoutInternal) = 1
QuantumOpticsBase.tensor(l::BravaisLattice, b::Basis) = Sample(l, b)
QuantumOpticsBase.tensor(b::Basis, l::BravaisLattice) = Sample(l, b)
QuantumOpticsBase.tensor(b::Basis, s::SampleWithoutInternal) = Sample(s.latt, s.boundaries, b)
QuantumOpticsBase.tensor(s::SampleWithoutInternal, b::Basis) = Sample(s.latt, s.boundaries, b)

@enum ParticleStatistics begin
    OneParticle = 0
    FermiDirac = 1
    BoseEinstein = -1
end

abstract type System{SampleT} end
Base.length(system::System) = length(system.sample)
sample(sys::System) = sys.sample
default_bonds(sys::System) = default_bonds(sys.sample)

struct FilledZones{SampleT} <: System{SampleT}
    sample::SampleT
    chempotential::Float64
    statistics::ParticleStatistics
    function FilledZones(sample::SampleT, μ; statistics=FermiDirac) where SampleT
        statistics == OneParticle && μ != 0 && error("OneParticle statistics cannot be present in systems with non-zero chemical potential")
        new{SampleT}(sample, μ, statistics)
    end
end

struct Particles{SampleT} <: System{SampleT}
    sample::SampleT
    nparticles::Int
    statistics::ParticleStatistics
    function Particles(sample::SampleT, nparticles; statistics=FermiDirac) where SampleT
        statistics == OneParticle && nparticles> 1 && error("OneParticle statistics not supported for multi-particle systems")
        new{SampleT}(sample, nparticles, statistics)
    end
end

function System(sample::Sample; μ = nothing, N = nothing, statistics=FermiDirac)
    if μ !== nothing && N === nothing
        return FilledZones(sample, μ, statistics=statistics)
    elseif N !== nothing && μ === nothing
        return Particles(sample, N, statistics=statistics)
    else
        error("Please define the chemical potential `μ` or the particle number `N` (but not both)")
    end
end
System(args...; μ = nothing, N = nothing, statistics=FermiDirac, kw...) =
    System(Sample(args...; kw...), μ=μ, N=N, statistics=statistics)

lattice(sys::System) = lattice(sys.sample)
Base.zero(sys::System) = zero(systembasis(sys))

function occupations(ps::Particles)
    if ps.statistics == FermiDirac
        fermionstates(length(ps.sample), ps.sample.nparticles)
    elseif ps.statistics == BoseEinstein
        bosonstates(length(ps.sample), ps.sample.nparticles)
    end
end

systembasis(sys::FilledZones) = basis(sys.sample)
systembasis(sys::Particles) = ManyBodyBasis(basis(sys.sample), occupations(sys))

function QuantumOpticsBase.manybodyoperator(ps::Particles, op::AbstractOperator)
    check_samebases(systembasis(ps), basis(op))
    return manybodyoperator(ManyBodyBasis(bas, occupations(ps)), oper)
end
function QuantumOpticsBase.manybodyoperator(ps::FilledZones, op::AbstractOperator)
    check_samebases(systembasis(ps), basis(op))
    return op
end
QuantumOpticsBase.manybodyoperator(sample, mat::AbstractMatrix) =
    manybodyoperator(sample, Operator(basis(sample), mat))

shift_site(sample::Sample, site) = shift_site(sample.boundaries, sample.latt, site)
shift_site(sys::System, site) = shift_site(sys.sample, site)

sample(lb::LatticeBasis) = Sample(lb.latt)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].latt, b.bases[1])
sample(b::Basis) = throw(MethodError(sample, (b,)))
sample(any) = sample(basis(any))
lattice(any) = lattice(sample(any))
internal_basis(sample::SampleWithInternal) = sample.internal
internal_basis(::SampleWithoutInternal) = throw(ArgumentError("Sample has no internal basis"))
internal_basis(any) = internal_basis(sample(any))
internal_length(sample::SampleWithInternal) = length(sample.internal)
internal_length(sample::SampleWithoutInternal) = 1
internal_length(any) = internal_length(sample(any))

QuantumOpticsBase.basis(sample::Sample) = sample.internal ⊗ LatticeBasis(sample.latt)
QuantumOpticsBase.basis(sample::SampleWithoutInternal) = LatticeBasis(sample.latt)
onebodybasis(sample::Sample) = basis(sample)
onebodybasis(sys::System) = onebodybasis(sys.sample)

macro accepts_system(fname, default_basis=nothing)
    esc(quote
        $fname(sample::Sample, args...; kw...) =
            $fname(FilledZones(sample, 0, statistics=FermiDirac), args...; kw...)
        $fname(l::BravaisLattice, bas::Basis, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, boundaries, bas), args...; kw...)
        $fname(l::BravaisLattice, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, boundaries, $default_basis), args...; kw...)
    end)
end

macro accepts_sample(fname, default_basis=nothing)
    esc(quote
        $fname(system::FilledZones, args...; kw...) =
            $fname(system.sample, args...; kw...)
        $fname(l::BravaisLattice, bas::Basis, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, boundaries, bas), args...; kw...)
        $fname(l::BravaisLattice, args...; boundaries=BoundaryConditions(), kw...) =
            $fname(Sample(l, boundaries, $default_basis), args...; kw...)
    end)
end
