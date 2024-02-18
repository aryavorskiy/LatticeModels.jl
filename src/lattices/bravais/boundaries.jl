import QuantumOpticsBase: Basis, AbstractOperator, basis, check_samebases

shift_site(js::SVector{N, Int}, lp::BravaisPointer{N}) where N =
    BravaisPointer(lp.latcoords + js, lp.basindex)
shift_site(js::SVector{N}, site::BravaisSite{N}) where N =
    BravaisSite(shift_site(js, bravaispointer(site)), site.unitcell)

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

function Base.show(io::IO, ::MIME"text/plain", bc::TwistedBoundary)
    print(io, bc.R, " → ")
    if bc.Θ % 2π ≈ 0
        print("periodic")
    else
        print("twist θ = ", trunc(bc.Θ, digits=2))
    end
end

struct FunctionBoundary{N, F<:Function} <: Boundary{N}
    condition::F
    R::SVector{N, Int}
end
FunctionBoundary(f::F, svec::AbstractArray) where F<:Function = FunctionBoundary{length(svec), F}(f, svec)
function shift_site(bc::FunctionBoundary{N}, i::Int, site::BravaisSite{N}) where N
    i == 0 && return 1., site
    factor = 1.
    for _ in 1:abs(i)
        if i > 0
            site = shift_site(-bc.R, site)
            factor *= bc.condition(site)
        else
            factor /= bc.condition(site)
            site = shift_site(bc.R, site)
        end
    end
    factor, site
end
Base.show(io::IO, ::MIME"text/plain", bc::FunctionBoundary) =
    print(io, bc.R, " → function '", bc.condition, "'")

to_boundary(p::Pair{<:AbstractVector, Bool}) =
    p[2] ? PeriodicBoundary(p[1]) : nothing
to_boundary(p::Pair{<:AbstractVector, <:Real}) =
    TwistedBoundary(p[1], p[2])
to_boundary(p::Pair{<:AbstractVector, <:Function}) =
    FunctionBoundary(p[2], p[1])

b_depth(::AbstractLattice) = 1
function route(bcs::BoundaryConditions{<:NTuple{M}}, l::AbstractLattice, lp::BravaisPointer{N}) where {M, N}
    for cind in cartesian_indices(b_depth(l), Val(M))
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
function shift_site(bcs::BoundaryConditions, l::AbstractLattice, site::BravaisSite)
    factor = 1.
    tup = route(bcs, l, bravaispointer(site))
    for i in 1:length(bcs.bcs)
        new_factor, site = shift_site(bcs.bcs[i], tup[i], site)
        factor *= new_factor
    end
    return factor, site
end
shift_site(::BoundaryConditions, ::AbstractLattice, ::NoSite) = 1., NoSite()

struct MagneticBoundaryConditions end
