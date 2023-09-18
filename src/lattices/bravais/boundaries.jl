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
TwistedBoundary(svec::AbstractArray, Θ::Real) = TwistedBoundary{length(svec)}(svec, Θ)
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
            site = shift_site(bc.R, l, site)
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
BoundaryConditions(::NTuple{M, <:Boundary} where M) =
    throw(ArgumentError("Dimension inconsistency. Check that cimension count for all boundaries is the same"))
BoundaryConditions(args...) =
    BoundaryConditions(Tuple(skipmissing(to_boundary.(args))))
to_boundary(p::Pair{<:AbstractVector, Bool}) =
    p[2] ? PeriodicBoundary(p[1]) : missing
to_boundary(p::Pair{<:AbstractVector, <:Real}) =
    TwistedBoundary(p[1], p[2])
to_boundary(p::Pair{<:AbstractVector, <:Function}) =
    FunctionBoundary(p[2], p[1])
to_boundary(b::Boundary) = b
to_boundary(b) = throw(ArgumentError("Could not convert `$b` to Boundary"))

@generated cartesian_indices(depth::Int, ::Val{M}) where M = quote
    CartesianIndex($((:(-depth) for _ in 1:M)...)):CartesianIndex($((:depth for _ in 1:M)...))
end
function route(bcs::BoundaryConditions{<:NTuple{M}}, l::AbstractLattice, lp::BravaisPointer{N}, depth=1) where {M, N}
    for cind in cartesian_indices(depth, Val(M))
        tup = Tuple(cind)
        tr_vec = @SVector zeros(Int, N)
        for i in 1:M
            tr_vec += tup[i] * bcs.bcs[i].R
        end
        new_lp = shift_site(-tr_vec, lp)
        new_lp in l && return tup
    end
    return Tuple(@SVector zeros(Int, M))
end
function shift_site(bcs::BoundaryConditions, l::AbstractLattice, site)
    factor = 1.
    tup = route(bcs, l, site.lp)
    for i in 1:length(bcs.bcs)
        new_factor, site = shift_site(bcs.bcs[i], tup[i], site)
        factor *= new_factor
    end
    return factor, site
end

struct MagneticBoundaryConditions end
