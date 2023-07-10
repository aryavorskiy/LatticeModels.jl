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

struct Sample{LT, BT, FT}
    latt::LT
    boundaries::BT
    field::FT
    nparticles::Int
end
Sample(latt, boundaries=BoundaryConditions(), field=NoField(); N=1) = Sample(latt, boundaries, field, N)
