abstract type AbstractSite{N} end

dims(::AbstractSite{N}) where N = N
Base.iterate(site::AbstractSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)

Base.show(io::IO, ::MIME"text/plain", site::AbstractSite{N}) where N =
    print(io, "Site of a ", N, "-dim lattice @ x = $(site.coords)")

struct NoSite <: AbstractSite{0} end
Base.show(io::IO, ::NoSite) = print(io, "LatticeModels.NoSite()")

# Site parameters
abstract type SiteProperty end
getsiteproperty(::AbstractSite, p::SiteProperty) = throw(ArgumentError("Site has no property `$p`"))

struct SitePropertyAlias{Symb} <: SiteProperty end
getsiteproperty(::AbstractSite, ::SitePropertyAlias{Symb}) where Symb =
    throw(ArgumentError("`$Symb` is not a valid site property"))
@inline getsiteproperty(site::AbstractSite, sym::Symbol) =
    getsiteproperty(site, SitePropertyAlias{sym}())

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

struct Coord <: SiteProperty axis::Int end
@inline function getsiteproperty(site::AbstractSite, c::Coord)
    @assert 1 ≤ c.axis ≤ dims(site)
    return getfield(site, :coords)[c.axis]
end
SitePropertyAlias{:x}() = Coord(1)
SitePropertyAlias{:y}() = Coord(2)
SitePropertyAlias{:z}() = Coord(3)
for i in 1:32
    @eval SitePropertyAlias{$(QuoteNode(Symbol("x$i")))}() = Coord($i)
end

abstract type AbstractLattice{SiteT} <: AbstractSet{SiteT} end
lattice(l::AbstractLattice) = l
dims(::AbstractLattice{<:AbstractSite{N}}) where {N} = N
dims(l) = dims(lattice(l))
Base.size(l::AbstractLattice) = (length(l),)
site_index(::AbstractLattice, ::NoSite) = nothing
Base.getindex(::AbstractLattice, ::Nothing) = NoSite()
Base.pairs(l::AbstractLattice) = Base.Iterators.Pairs(l, 1:length(l))

Base.summary(io::IO, l::AbstractLattice{<:AbstractSite{N}}) where N =
    print(io, length(l), "-site ", N, "-dim lattice")
function Base.show(io::IO, mime::MIME"text/plain", l::AbstractLattice)
    summary(io, l)
    if !get(io, :compact, false)
        print(io, ":")
        for i in 1:min(length(l), 10)
            print("\n  ")
            show(io, mime, l[i])
        end
        length(l) > 10 && print("\n  ⋮")
    end
end

# Set functions
Base.copy(l::AbstractLattice) = Base.copymutable(l)
Base.in(site::SiteT, l::AbstractLattice{SiteT}) where {SiteT} =
    site_index(l, site) !== nothing
Base.hasfastin(::Type{<:AbstractLattice}) = true
function Base.delete!(l::AbstractLattice{ST}, site::ST) where ST<:AbstractSite
    i = site_index(l, site)
    i !== nothing && Base.deleteat!(l, i)
    return l
end

# Indexing
Base.firstindex(::AbstractLattice) = 1
Base.lastindex(l::AbstractLattice) = length(l)
Base.eachindex(l::AbstractLattice) = firstindex(l):lastindex(l)
Base.checkbounds(::Type{Bool}, l::AbstractLattice, is) = all(1 .≤ is .≤ length(l))
Base.checkbounds(l::AbstractLattice, is) =
    !Base.checkbounds(Bool, l, is) && throw(BoundsError(l, is))
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
    sites(l1) != sites(l2) &&
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
    index::Int
    factor::Float64
end
function ResolvedSite(site::ST, index::Int) where {ST}
    ResolvedSite{ST}(site, index, 1.0)
end
function resolve_site(l::AbstractLattice, site::AbstractSite)
    index = site_index(l, site)
    index === nothing && return nothing
    ResolvedSite(site, index)
end

const SingleBond{LT<:AbstractSite} = Pair{LT, LT}
