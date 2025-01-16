"""
    AbstractSite{N}

An abstract type for a site of a `N`-dimensional lattice.

## Fields
- `coords`: A `SVector` of size `N` representing the spatial coordinates of the site.
    All subtypes are expected to have this field.
"""
abstract type AbstractSite{N} end

dims(::AbstractSite{N}) where N = N
Base.iterate(site::AbstractSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)

Base.summary(io::IO, site::AbstractSite{N}) where N =
    print(io, N, "-dim ", typeof(site))
function Base.show(io::IO, ::MIME"text/plain", site::AbstractSite{N}) where N
    if site in get(io, :SHOWN_SET, ())
        print(io, "Site at ", site.coords)
    else
        summary(io, site)
        print(io, " at ", site.coords)
    end
end

struct NoSite <: AbstractSite{0} end
Base.show(io::IO, ::MIME"text/plain", ::NoSite) = print(io, "LatticeModels.NoSite()")

"""
    SiteProperty

An abstract type for a property of a site.

This interface is used to define various properties of a site. They can be accessed using
`getsiteproperty`. This interface is used in following places:
- `lattice[...]` syntax to access sites with specific properties.
- `lattice_value[...]` syntax to access values defined on sites with specific properties.
- Functions to generate `LatticeValue`s and operators for specific properties.

## Examples
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(3, 3);

julia> l[x = 1, y = 2]          # Get site with x = 1 and y = 2
2-dim Bravais lattice site in 2D space at [1.0, 2.0]

julia> l[x = 1]                 # Get sublattice with x = 1
3-site 2-dim Bravais lattice in 2D space
Unit cell:
  Basis site coordinates:
    ┌      ┐
    │ 0.000│
    │ 0.000│
    └      ┘
  Translation vectors:
    ┌      ┐ ┌      ┐
    │ 1.000│ │ 0.000│
    │ 0.000│ │ 1.000│
    └      ┘ └      ┘
Lattice type: SquareLattice{2}
Default translations:
  :axis1 → Bravais[3, 0]
  :axis2 → Bravais[0, 3]
Nearest neighbor hoppings:
  1.00000 =>
    Bravais[1, 0]
    Bravais[0, 1]
  1.41421 =>
    Bravais[1, -1]
    Bravais[1, 1]
  2.00000 =>
    Bravais[2, 0]
    Bravais[0, 2]
Boundary conditions: none

julia> l[x = 1, y = 2, z = 3]   # No site with defined z property on a 2D lattice
ERROR: ArgumentError: Invalid axis index 3 of a 2-dim site
[...]
```
"""
abstract type SiteProperty end
getsiteproperty(::AbstractSite, p::SiteProperty) = throw(ArgumentError("Site has no property `$p`"))

struct SitePropertyAlias{Symb} <: SiteProperty end
getsiteproperty(::AbstractSite, ::SitePropertyAlias{Symb}) where Symb =
    throw(ArgumentError("`$Symb` is not a valid site property"))
@inline getsiteproperty(site::AbstractSite, sym::Symbol) =
    getsiteproperty(site, SitePropertyAlias{sym}())

@static if VERSION ≥ v"1.8"
    @inline function Base.getproperty(site::AbstractSite, sym::Symbol)
        if sym in fieldnames(typeof(site))
            return getfield(site, sym)
        end
        prop = SitePropertyAlias{sym}()
        if !(prop isa SitePropertyAlias)
            # If the constructor was overridden, then get site property
            return getsiteproperty(site, prop)
        end
        return getfield(site, sym)
    end
end

struct Coord <: SiteProperty axis::Int end
@inline function getsiteproperty(site::AbstractSite, c::Coord)
    1 ≤ c.axis ≤ dims(site) ||
        throw(ArgumentError("Invalid axis index $(c.axis) of a $(dims(site))-dim site"))
    return getfield(site, :coords)[c.axis]
end
SitePropertyAlias{:x}() = Coord(1)
SitePropertyAlias{:y}() = Coord(2)
SitePropertyAlias{:z}() = Coord(3)
for i in 1:32
    @eval SitePropertyAlias{$(QuoteNode(Symbol("x$i")))}() = Coord($i)
end

"""
    AbstractLattice{SiteT}

An abstract type for a lattice of `SiteT` sites.

## Methods for subtypes to implement
- `length(l::AbstractLattice)`: Return the number of sites in the lattice.
- `site_index(l::AbstractLattice, site::SiteT)`: Return the index of the site in the lattice.
- `getindex(l::AbstractLattice, i::Int)`: Return the site with the given index.
- `getindex(l::AbstractLattice, is::AbstractVector{Int})`: Return an `AbstractLattice` with the sites at the given indices.

## Optional methods for mutable lattices
- `emptymutable(l::AbstractLattice, ::Type{SiteT})`: Return an empty mutable instance of lattice.
- `copymutable(l::AbstractLattice)`: Return a mutable copy of the lattice.
- `push!(l::AbstractLattice, site::SiteT)`: Add a site to the lattice.
- `deleteat!(l::AbstractLattice, is::AbstractVector{Int})`: Remove the sites with the given indices from the lattice.
"""
abstract type AbstractLattice{SiteT} <: AbstractSet{SiteT} end

"""
    lattice(any)

Return the lattice of the given object (an operator, `LatticeValue`, ...)
"""
lattice(l::AbstractLattice) = l
dims(::AbstractLattice{<:AbstractSite{N}}) where {N} = N
dims(l) = dims(lattice(l))
Base.size(l::AbstractLattice) = (length(l),)

"""
    site_index(lat, site[, range])

Return the index of the `site` in the lattice `lat`. If `range` is given, only search in the
given range. Return `nothing` if the site is not found.
"""
site_index(::AbstractLattice, ::NoSite, range=nothing) = nothing
site_index(l::AbstractLattice, site::AbstractSite) = site_index(l, site, eachindex(l))
Base.getindex(::AbstractLattice, ::Nothing) = NoSite()
Base.pairs(l::AbstractLattice) = Base.Iterators.Pairs(l, 1:length(l))

latticename(any) = string(typeof(any))
Base.summary(io::IO, l::AbstractLattice{<:AbstractSite{N}}) where N =
    print(io, length(l), "-site ", latticename(l), " in ", N, "D space")
function Base.show(io::IO, mime::MIME"text/plain", l::AbstractLattice)
    summary(io, l)
    if !requires_compact(io) && length(l) > 0
        io = IOContext(io, :compact => true, :SHOWN_SET => l)
        print(io, ":")
        maxlen = get(io, :maxlines, 10)
        for i in 1:min(length(l), maxlen)
            print(io, "\n  ")
            if i == maxlen < length(l)
                print(io, "  ⋮")
            else
                show(io, mime, l[i])
            end
        end
    end
end

# Set functions
Base.copy(l::AbstractLattice) = Base.copymutable(l)
Base.in(site::SiteT, l::AbstractLattice{SiteT}) where {SiteT} =
    site_index(l, site) !== nothing
function Base.delete!(l::AbstractLattice{ST}, site::ST) where ST<:AbstractSite
    i = site_index(l, site)
    i !== nothing && Base.deleteat!(l, i)
    return l
end

# iteration (assume that the lattice has fast indexing - usually it does)
function Base.iterate(l::AbstractLattice, state = (1, length(l)))
    i, len = state
    return i > len ? nothing : (l[i], (i+1, len))
end

# Indexing
Base.firstindex(::AbstractLattice) = 1
Base.lastindex(l::AbstractLattice) = length(l)
Base.eachindex(l::AbstractLattice) = firstindex(l):lastindex(l)
Base.checkbounds(::Type{Bool}, l::AbstractLattice, is::Union{Int,AbstractVector{Int}}) =
    all(1 .≤ is .≤ length(l))
Base.checkbounds(l::AbstractLattice, is) =
    !checkbounds(Bool, l, is) && throw(BoundsError(l, is))
Base.push!(l::AbstractLattice{SiteT}, ::SiteT) where SiteT = error("Define `push!` for $l")
Base.push!(l::AbstractLattice{SiteT}, x::Any) where SiteT = push!(l, convert(SiteT, x))
function Base.push!(l::AbstractLattice, xs::Any...)
    for x in xs
        push!(l, x)
    end
    return l
end
Base.pop!(l::AbstractLattice) = deleteat!(l, lastindex(l))
Base.popfirst!(l::AbstractLattice) = deleteat!(l, firstindex(l))

function Base.filter!(f::Function, l::AbstractLattice)
    is = Int[]
    for (i, site) in enumerate(l)
        f(site) || push!(is, i)
    end
    deleteat!(l, is)
    return l
end

function check_param_pairs(pairs::Tuple{Vararg{Pair}})
    for (prop, _) in pairs
        if prop isa Symbol
            error("`:$prop => ...` notation is purposely disallowed. Use `$prop = ...`")
        end
    end
end
kw_to_param_pairs(ntup::NamedTuple{T}) where T =
    tuple(SitePropertyAlias{first(T)}() => first(ntup),
        kw_to_param_pairs(Base.structdiff(ntup, NamedTuple{(first(T),)}))...)
kw_to_param_pairs(::NamedTuple{()}) = ()
kw_to_param_pairs(ps::Base.Iterators.Pairs) = kw_to_param_pairs(NamedTuple(ps))

function to_param_pairs(arg...; kw...)
    check_param_pairs(arg)
    return tuple(arg..., kw_to_param_pairs(kw)...)
end

match_param_pairs(site, ::Tuple{}) = true
function match_param_pairs(site, pairs::Tuple{Vararg{Pair}})
    param, val = first(pairs)
    return getsiteproperty(site, param) in val && match_param_pairs(site, Base.tail(pairs))
end

function pairs_to_index(l::AbstractLattice, all_pairs::Tuple{Vararg{Pair}})
    ind = 0
    nfound = 0
    for (i, site) in enumerate(l)
        if match_param_pairs(site, all_pairs)
            ind = i
            nfound += 1
        end
    end
    nfound > 1 && throw(ArgumentError("More than one site satisfies parameter conditions"))
    nfound < 1 && return nothing
    return ind
end
function pairs_to_indices(l::AbstractLattice, all_pairs::Tuple{Vararg{Pair}})
    inds = Int[]
    for (i, site) in enumerate(l)
        if match_param_pairs(site, all_pairs)
            push!(inds, i)
        end
    end
    return inds
end

function collect_coords(l::AbstractLattice)
    d = dims(l)
    pts = zeros(d, length(l))
    for (i, site) in enumerate(l)
        pts[:, i] = site.coords
    end
    pts
end

"""
    IncompatibleLattices([header, ]lat1, lat2)

An exception thrown when two lattices are incompatible.
"""
struct IncompatibleLattices <: Exception
    header::String
    l1::AbstractLattice
    l2::AbstractLattice
    IncompatibleLattices(header, l1, l2) = new(header, lattice(l1), lattice(l2))
end

function Base.showerror(io::IO, ex::IncompatibleLattices)
    print(io, "$(ex.header).\nGot following:\n  #1: ")
    io = IOContext(io, :compact => true)
    show(io, "text/plain", ex.l1)
    print(io, "\n  #2: ")
    show(io, "text/plain", ex.l2)
end

"""
Checks if `l1` and `l2` objects are defined on the same lattice. Throws an error if not.
"""
function check_samelattice(l1, l2)
    lattice(l1) !== lattice(l2) && lattice(l1) != lattice(l2) &&
        throw(IncompatibleLattices("Matching lattices expected", l1, l2))
end

"""
Checks if `l1` and `l2` objects are defined on the same sites. Throws an error if not.
"""
function check_samesites(l1, l2)
    stripmeta(lattice(l1)) != stripmeta(lattice(l2)) &&
        throw(IncompatibleLattices("Matching sets of sites expected", l1, l2))
end

"""
Checks if `l1` is sublattice of `l2`. Throws an error if not.
"""
function check_issublattice(l1::AbstractLattice, l2::AbstractLattice)
    !issubset(l1, l2) &&
        throw(IncompatibleLattices("#1 is expected to be sublattice of #2", l1, l2))
end

struct ResolvedSite{ST}
    site::ST
    old_site::ST
    index::Int
    factor::ComplexF64
    function ResolvedSite(site::ST, old_site::ST, index::Int, factor=1) where ST<:AbstractSite
        new{ST}(site, old_site, index, ComplexF64(factor))
    end
end
ResolvedSite(site::AbstractSite, index::Int) = ResolvedSite(site, site, index)
@inline function resolve_site_default(l::AbstractLattice, site::AbstractSite)
    index = site_index(l, site)
    index === nothing && return nothing
    ResolvedSite(site, index)
end
resolve_site(l::AbstractLattice, site::AbstractSite) = resolve_site_default(l, site)
resolve_site(::AbstractLattice, rs::ResolvedSite) = rs
resolve_site(l::AbstractLattice, i::Int) =
    i in eachindex(l) ? ResolvedSite(l[i], i) : nothing

struct LatticeWithMetadata{LT,MetaT,SiteT} <: AbstractLattice{SiteT}
    lat::LT
    metadata::MetaT
    function LatticeWithMetadata(lat::LT, meta::MetaT) where
            {SiteT,LT<:AbstractLattice{SiteT},MetaT<:NamedTuple}
        new{LT,MetaT,SiteT}(lat, meta)
    end
end
LatticeWithMetadata(lw::LatticeWithMetadata, newmeta::NamedTuple) =
    LatticeWithMetadata(lw.lat, merge(lw.metadata, newmeta))
allmeta(::AbstractLattice) = NamedTuple()
allmeta(lw::LatticeWithMetadata) = lw.metadata

Base.:(==)(lw1::LatticeWithMetadata, lw2::LatticeWithMetadata) = lw1.lat == lw2.lat
Base.length(lw::LatticeWithMetadata) = length(lw.lat)
Base.getindex(::LatticeWithMetadata, ::Nothing) = NoSite()
Base.getindex(lw::LatticeWithMetadata, i::Int) = lw.lat[i]
Base.getindex(lw::LatticeWithMetadata, is::AbstractVector{Int}) =
    LatticeWithMetadata(lw.lat[is], lw.metadata)
site_index(::LatticeWithMetadata, ::NoSite) = nothing

Base.emptymutable(l::LatticeWithMetadata, ::Type{T}) where {T<:AbstractSite} =
    LatticeWithMetadata(Base.emptymutable(l.lat, T), l.metadata)
Base.copymutable(lw::LatticeWithMetadata) =
    LatticeWithMetadata(Base.copymutable(lw.lat), lw.metadata)
Base.deleteat!(lw::LatticeWithMetadata, is) = (deleteat!(lw.lat, is); return lw)
function Base.push!(lw::LatticeWithMetadata{<:AbstractLattice{SiteT}}, site::SiteT) where SiteT
    push!(lw.lat, site)
    return lw
end

Base.show(io::IO, ::Type{<:LatticeWithMetadata{LT, MetaT}}) where {LT, MetaT} =
    print(io, "LatticeWithMetadata{", LT, ", params", fieldnames(MetaT), "}")
latticename(lw::LatticeWithMetadata) =
    hasmeta(lw, :latticetype) ? latticename(lw.latticetype) : latticename(lw.lat)
function Base.show(io::IO, mime::MIME"text/plain", lw::LatticeWithMetadata)
    io = IOContext(io, :maxlines=>4)
    requires_compact(io) && return summary(io, lw)
    show(io, mime, lw.lat)
    for v in values(lw.metadata)
        println(io)
        show(io, mime, v)
    end
end

function Base.getproperty(lw::LatticeWithMetadata, sym::Symbol)
    if sym in fieldnames(typeof(lw))
        return getfield(lw, sym)
    end
    meta = getfield(lw, :metadata)
    if haskey(meta, sym)
        return meta[sym]
    end
    return getproperty(getfield(lw, :lat), sym)
end

Base.propertynames(lw::LatticeWithMetadata) =
    tuple(propertynames(lw.lat)..., fieldnames(typeof(lw))..., keys(lw.metadata)...)

LatticeWithMetadata(l::AbstractLattice; kw...) = LatticeWithMetadata(l, NamedTuple(kw))

hasmeta(::AbstractLattice, ::Symbol) = false
hasmeta(l::LatticeWithMetadata, meta::Symbol) = haskey(l.metadata, meta)
getmeta(::AbstractLattice, ::Symbol, default=nothing) = default
getmeta(l::LatticeWithMetadata, meta::Symbol, default=nothing) = get(l.metadata, meta, default)
setmeta(l::AbstractLattice, param::Symbol, val) = LatticeWithMetadata(l, NamedTuple{(param,)}((val,)))
delmeta(l::LatticeWithMetadata, param::Symbol) = LatticeWithMetadata(l.lat, Base.structdiff(l.metadata, NamedTuple{(param,)}))

stripmeta(l::AbstractLattice) = l
stripmeta(l::LatticeWithMetadata) = l.lat

pushmeta(l::AbstractLattice, param::Symbol, val) =
    LatticeWithMetadata(stripmeta(l), merge(NamedTuple{(param,)}((val,)), allmeta(l)))

const MaybeWithMetadata{LT} = Union{LT, LatticeWithMetadata{<:LT}}
