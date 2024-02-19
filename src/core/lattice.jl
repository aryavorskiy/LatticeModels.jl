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

function Base.show(io::IO, ::MIME"text/plain", site::AbstractSite{N}) where N
    if site in get(io, :SHOWN_SET, ())
        print(io, "$(site.coords)")
    else
        print(io, N, "-dim ", typeof(site), " @ x = $(site.coords)")
    end
end

struct NoSite <: AbstractSite{0} end
Base.show(io::IO, ::NoSite) = print(io, "LatticeModels.NoSite()")

"""
    SiteProperty

An abstract type for a property of a site.

This interface is used to define various properties of a site. They can be accessed using
`getsiteproperty`. This interface is used in following places:
- `lattice[...]` syntax to access sites with specific properties.
- `lattice_value[...]` syntax to access values defined on sites with specific properties.
- `siteproperty_value` function to generate `LatticeValue` objects for specific properties.
- `siteproperty_operator` function to generate operators for specific properties.

## Examples
```jldoctest
julia> l = SquareLattice(3, 3)

julia> l[x = 1, y = 2]          # Get site with x = 1 and y = 2
Site of a 2-dim lattice @ [1.0, 2.0]

julia> l[x = 1]                 # Get sublattice with x = 1
3-site 2-dim SquareLattice

julia> l[x = 1, y = 2, z = 3]   # No site with defined z property on a 2D lattice
ERROR: ArgumentError: Invalid axis index 3 of a 2-dim site
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
lattice(l::AbstractLattice) = l
dims(::AbstractLattice{<:AbstractSite{N}}) where {N} = N
dims(l) = dims(lattice(l))
Base.size(l::AbstractLattice) = (length(l),)
site_index(::AbstractLattice, ::NoSite) = nothing
Base.getindex(::AbstractLattice, ::Nothing) = NoSite()
Base.pairs(l::AbstractLattice) = Base.Iterators.Pairs(l, 1:length(l))

Base.summary(io::IO, l::AbstractLattice{<:AbstractSite{N}}) where N =
    print(io, length(l), "-site ", N, "-dim ", typeof(l))
function Base.show(io::IO, mime::MIME"text/plain", l::AbstractLattice)
    summary(io, l)
    if !requires_compact(io) && length(l) > 0
        io = IOContext(io, :compact => true)
        print(io, ":")
        for i in 1:min(length(l), 10)
            print(io, "\n  ")
            show(io, mime, l[i])
        end
        length(l) > 10 && print(io, "\n  ⋮")
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
Base.pop!(l::AbstractLattice) = Base.deleteat!(l, lastindex(l))
Base.popfirst!(l::AbstractLattice) = Base.deleteat!(l, firstindex(l))

sym_to_param_pair(p::Pair{Symbol}) = SitePropertyAlias{p[1]}() => p[2]
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

function Base.getindex(l::AbstractLattice, pairs::Pair...; kw...)
    inds = pairs_to_indices(l, to_param_pairs(pairs...; kw...))
    length(inds) == 1 ? l[only(inds)] : l[inds]
end
Base.getindex(l::AbstractLattice, ::typeof(!), pairs::Pair...; kw...) =
    l[pairs_to_index(l, to_param_pairs(pairs...; kw...))]

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
IncompatibleLattices(l1, l2) = IncompatibleLattices("Incompatible lattices", l1, l2)

Base.showerror(io::IO, ex::IncompatibleLattices) = print(io,
"""$(ex.header).\nGot following:
        #1: $(repr("text/plain", ex.l1))
        #2: $(repr("text/plain", ex.l2))""")

"""
Checks if `l1` and `l2` objects are defined on the same lattice. Throws an error if not.
"""
function check_samelattice(l1, l2)
    lattice(l1) != lattice(l2) &&
        throw(IncompatibleLattices("Matching lattices expected", l1, l2))
end

"""
Checks if `l1` and `l2` objects are defined on the same sites. Throws an error if not.
"""
function check_samesites(l1, l2)
    stripparams(lattice(l1)) != stripparams(lattice(l2)) &&
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
    function ResolvedSite(site::ST, old_site::ST, index::Int, factor) where ST<:AbstractSite
        new{ST}(site, old_site, index, ComplexF64(factor))
    end
end
ResolvedSite(site::AbstractSite, old_site::AbstractSite, index::Int) =
    ResolvedSite(site, old_site, index, 1)
ResolvedSite(site::AbstractSite, index::Int) = ResolvedSite(site, site, index, 1)
function resolve_site(l::AbstractLattice, site::AbstractSite)
    index = site_index(l, site)
    index === nothing && return nothing
    ResolvedSite(site, index)
end

struct LatticeWithParams{LT,ParamsT,SiteT} <: AbstractLattice{SiteT}
    lat::LT
    params::ParamsT
    function LatticeWithParams(latt::LT, params::ParamsT) where
            {SiteT,LT<:AbstractLattice{SiteT},ParamsT<:NamedTuple}
        new{LT,ParamsT,SiteT}(latt, params)
    end
end
LatticeWithParams(lw::LatticeWithParams, params::NamedTuple) =
    LatticeWithParams(lw.lat, merge(lw.params, params))

Base.length(lw::LatticeWithParams) = length(lw.lat)
Base.getindex(::LatticeWithParams, ::Nothing) = NoSite()
Base.getindex(lw::LatticeWithParams, i::Int) = lw.lat[i]
Base.getindex(lw::LatticeWithParams, is::AbstractVector{Int}) = LatticeWithParams(lw.lat[is], lw.params)
site_index(lw::LatticeWithParams, site::AbstractSite) = site_index(lw.lat, site)
site_index(::LatticeWithParams, ::NoSite) = nothing

Base.emptymutable(l::LatticeWithParams, ::Type{T}) where {T <: AbstractSite} =
    LatticeWithParams(Base.emptymutable(l.lat, T), l.params)
Base.copymutable(lw::LatticeWithParams) = LatticeWithParams(copymutable(lw.lat), lw.params)
Base.deleteat!(lw::LatticeWithParams, is) = (deleteat!(lw.lat, is); return lw)
Base.push!(lw::LatticeWithParams, site) = (push!(lw.lat, site); return lw)

Base.summary(io::IO, lw::LatticeWithParams) = summary(io, lw.lat)
function Base.show(io::IO, mime::MIME"text/plain", lw::LatticeWithParams)
    if requires_compact(io)
        print(io, "Wrapped ")
        return show(io, mime, lw.lat)
    end
    show(io, mime, lw.lat)
    io = IOContext(io, :showtitle => false)
    for v in values(lw.params)
        println(io)
        show(io, mime, v)
    end
end

function Base.getproperty(lw::LatticeWithParams, sym::Symbol)
    if sym in fieldnames(typeof(lw))
        return getfield(lw, sym)
    end
    params = getfield(lw, :params)
    if haskey(params, sym)
        return params[sym]
    end
    return getproperty(getfield(lw, :lat), sym)
end

Base.propertynames(lw::LatticeWithParams) =
    tuple(propertynames(lw.lat)..., fieldnames(typeof(lw))..., keys(lw.params)...)

const Lattice = LatticeWithParams
Lattice(l::AbstractLattice; kw...) = LatticeWithParams(l, NamedTuple(kw))

hasparam(::AbstractLattice, ::Symbol) = false
hasparam(l::LatticeWithParams, param::Symbol) = haskey(l.params, param)
getparam(::AbstractLattice, ::Symbol, default=nothong) = default
getparam(l::LatticeWithParams, param::Symbol, default=nothing) = get(l.params, param, default)
setparam(l::AbstractLattice, param::Symbol, val) = Lattice(l, NamedTuple{(param,)}((val,)))

stripparams(l::AbstractLattice) = l
stripparams(l::LatticeWithParams) = l.lat

const OnSites{LT} = Union{LT, LatticeWithParams{<:LT}}
