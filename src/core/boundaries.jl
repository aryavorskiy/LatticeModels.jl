abstract type Boundary{TranslationT} end

to_boundary(any) = throw(ArgumentError("Could not convert `$any` to a `Boundary`"))
to_boundary(::Nothing) = nothing
to_boundary(b::Boundary) = b

"""
    BoundaryConditions

A collection of boundary conditions for a lattice.

## Fields
- `bcs`: A tuple of boundary conditions.
- `depth`: The upper limit of the depth of the boundary conditions (used for routing).
"""
struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    depth::Int
    function BoundaryConditions(bcs::CondsTuple; depth=2) where CondsTuple<:Tuple{Vararg{Boundary}}
        new{CondsTuple}(bcs, depth)
    end
end

_skipnothing(args::Tuple) = _skipnothing((), args)
_skipnothing(pre_args::Tuple, post_args::Tuple) =
    first(post_args) === nothing ?
    _skipnothing(pre_args, Base.tail(post_args)) :
    _skipnothing((pre_args..., first(post_args)), Base.tail(post_args))
_skipnothing(pre_args::Tuple, ::Tuple{}) = pre_args

BoundaryConditions(args...; depth=2) =
    BoundaryConditions(_skipnothing(to_boundary.(args)); depth=depth)

function to_boundaries(arg)
    if arg isa Tuple
        return BoundaryConditions(arg...)
    elseif arg isa Boundary
        return BoundaryConditions(arg)
    elseif arg isa BoundaryConditions
        return arg
    else
        return BoundaryConditions(to_boundary(arg))
    end
end

getboundaries(l::AbstractLattice) = getparam(l, :boundaries, BoundaryConditions())
setboundaries(l::AbstractLattice, bcs::BoundaryConditions) = setparam(l, :boundaries, bcs)
setboundaries(l::AbstractLattice, bcs) = setboundaries(l, to_boundaries(bcs))

Base.getindex(bcs::BoundaryConditions, i::Int) = bcs.bcs[i]

function Base.show(io::IO, mime::MIME"text/plain", bcs::BoundaryConditions)
    print(io, "Boundary conditions:")
    for bc in bcs.bcs
        println()
        show(io, mime, bc)
    end
end

function resolve_site(l::LatticeWithParams, site::AbstractSite)
    factor, site = shift_site(getboundaries(l), l, site)
    i = site_index(l, site)
    i === nothing && return nothing
    return ResolvedSite(site, i, factor)
end

cartesian_indices(b::BoundaryConditions{NTuple{M}}) where M =
    cartesian_indices(b.b_depth, Val{M}())
cartesian_indices(l) = cartesian_indices(getboundaries(l))
