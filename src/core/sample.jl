abstract type Boundary end

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
    site = shift_site(l, site, -dv * offset * lspan)
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

struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    function BoundaryConditions(bcs::CondsTuple) where CondsTuple<:NTuple{N, <:Boundary} where N
        @assert allunique(bc.i for bc in bcs)
        new{CondsTuple}(bcs)
    end
end
_extract_boundary_conditions(b::Boundary) = b
function _extract_boundary_conditions(pb::Pair{Int, Bool})
    !pb.second && error("")
    PeriodicBoundary(pb.first)
end
_extract_boundary_conditions(pb::Pair{Int, <:Real}) = TwistedBoundary(pb...)
BoundaryConditions(args...) = BoundaryConditions(_extract_boundary_conditions.(args))
function shift_site(bcs::BoundaryConditions, l::Lattice, site::LatticeSite)
    factor = 1.
    for bc in bcs.bcs
        new_factor, site = shift_site(bc, l, site)
        factor *= new_factor
    end
    factor, site
end

@enum ParticleStatistics begin
    one_particle
    fermi
    bose
end

struct Sample{LT, BasisT}
    latt::LT
    internal::BasisT
    nparticles::Int
    statistics::ParticleStatistics
end

function Sample(latt::LT, internal::BT=GenericBasis(1);
        N::Int=1, statistics::ParticleStatistics=one_particle) where {LT, BT}
    N ≤ 0 && error("Positive particle count expected")
    N == 1 && return Sample(latt, internal, 1, one_particle)
    statistics == one_particle && error("One-particle statistics invalid for multi-particle systems")
    Sample{LT, BT}(latt, internal, N, statistics)
end
Base.length(sample::Sample) = length(sample.latt) * length(sample.internal)
lattice(sample::Sample) = sample.latt
internal_one(sample::Sample) =
    sample.statistics == one_particle ? 1 : one(sample.internal)
