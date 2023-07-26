import QuantumOpticsBase: Basis, AbstractOperator
abstract type Boundary end
abstract type AbstractBoundaryConditions end

struct TwistedBoundary <: Boundary
    i::Int
    Θ::Float64
end
PeriodicBoundary(i) = TwistedBoundary(i, 0)
function shift_site(bc::TwistedBoundary, l::Lattice, site::LatticeSite{N}) where N
    ret = 1., site
    bc.i > dims(l) && return ret
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return ret
    dv = one_hot(bc.i, Val(N))
    site = shift_site(-dv * offset * lspan, l, site)
    exp(im * bc.Θ * offset), site
end

struct FunctionBoundary{F<:Function} <: Boundary
    condition::F
    i::Int
end

shift_site(js::SVector{N}, l::Lattice, site::LatticeSite{N}) where N =
    LatticeSite(site.unit_cell + js, site.basis_index, site.coords + bravais(l).translation_vectors * js)

function shift_site(bc::FunctionBoundary, l::Lattice, site::LatticeSite{N}) where N
    ret = 1., site
    bc.i > dims(l) && return ret
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return ret
    dv = one_hot(bc.i, Val(N))
    factor = 1.
    for _ in 1:abs(offset)
        if offset > 0
            site = shift_site(-dv * lspan, l, site)
            factor *= bc.condition(ls)
        else
            factor *= bc.condition(ls)'
            site = shift_site(dv * lspan, l, site)
        end
    end
    factor, site
end

struct BoundaryConditions{CondsTuple} <: AbstractBoundaryConditions
    bcs::CondsTuple
    function BoundaryConditions(bcs::CondsTuple) where CondsTuple<:NTuple{N, <:Boundary} where N
        @assert allunique(bc.i for bc in bcs)
        new{CondsTuple}(bcs)
    end
end

struct MagneticBoundaryConditions <: AbstractBoundaryConditions end
_extract_boundary_conditions(b::Boundary) = b
function _extract_boundary_conditions(pb::Pair{Int, Bool})
    !pb.second && return missing
    PeriodicBoundary(pb.first)
end
_extract_boundary_conditions(pb::Pair{Int, <:Real}) = TwistedBoundary(pb...)
BoundaryConditions(args...) = BoundaryConditions(Tuple(skipmissing(_extract_boundary_conditions.(args))))
function shift_site(bcs::BoundaryConditions, l::Lattice, site::LatticeSite)
    factor = 1.
    for bc in bcs.bcs
        new_factor, site = shift_site(bc, l, site)
        factor *= new_factor
    end
    factor, site
end

@enum ParticleStatistics begin
    OneParticle = 0
    FermiDirac = 1
    BoseEinstein = -1
end

struct Sample{AdjMatT, LT, BasisT, BoundaryT}
    adjacency_matrix::AdjMatT
    latt::LT
    boundaries::BoundaryT
    internal::BasisT
    nparticles::Int
    statistics::ParticleStatistics
end
const LatticeSample{AdjMatT, LT, BoundaryT} = Sample{AdjMatT, LT, Nothing, BoundaryT}

function Sample(adjacency_matrix, latt::LT, internal::BT=nothing; N::Int=1,
        statistics::ParticleStatistics=OneParticle, boundaries=BoundaryConditions()) where {LT<:Lattice, BT<:Nullable{Basis}}
    N ≤ 0 && error("Positive particle count expected")
    N != 1 && statistics == OneParticle && error("One-particle statistics invalid for multi-particle systems")
    return Sample(adjacency_matrix, latt, boundaries, internal, N, statistics)
end
Sample(latt::Lattice, internal=nothing; kw...) =
    Sample(nothing, latt, internal; kw...)

Base.length(sample::Sample) = length(sample.latt) * length(sample.internal)
Base.length(sample::LatticeSample) = length(sample.latt)
lattice(sample::Sample) = sample.latt
default_bonds(sample::Sample) = default_bonds(lattice(sample))
internal_one(sample::Sample) = one(sample.internal)
internal_one(sample::LatticeSample) = 1
ismanybody(sample::Sample) = sample.nparticles != 1

function occupations(sample::Sample)
    sample.nparticles == 1 && throw(ArgumentError("Cannot generate occupations for one-particle sample"))
    if sample.statistics == FermiDirac
        fermionstates(length(sample), sample.nparticles)
    elseif sample.statistics == BoseEinstein
        bosonstates(length(sample), sample.nparticles)
    end
end

samplebasis(sample::Sample) = sample.nparticles == 1 ? onebodybasis(sample) :
        ManyBodyBasis(onebodybasis(sample), occupations(sample))

QuantumOpticsBase.:(⊗)(l::Lattice, b::Basis) = Sample(l, b)
QuantumOpticsBase.:(⊗)(b::Basis, l::Lattice) = Sample(l, b)
Base.zero(sample::Sample) = zero(samplebasis(sample))

function QuantumOpticsBase.manybodyoperator(sample::Sample, op::AbstractOperator)
    if sample.statistics == OneParticle
        return op
    else
        return manybodyoperator(ManyBodyBasis(bas, occupations(sample)), oper)
    end
end
QuantumOpticsBase.manybodyoperator(sample::Sample, mat::AbstractMatrix) =
    manybodyoperator(sample, Operator(onebodybasis(sample), mat))
