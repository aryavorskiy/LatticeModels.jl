"""
    GenericSite{N}

A generic site in an N-dimensional lattice.
"""
struct GenericSite{N} <: AbstractSite{N}
    coords::SVector{N,Float64}
    function GenericSite(coord::SVector{N,<:Real}) where N
        new{N}(convert(SVector{N,Float64}, coord))
    end
end
GenericSite(args...) = GenericSite(SVector(args))
Base.:(==)(s1::GenericSite, s2::GenericSite) = s1.coords == s2.coords
Base.isless(s1::GenericSite, s2::GenericSite) = isless(s1.coords, s2.coords)

Base.convert(::Type{GenericSite}, vec::AbstractVector{<:Real}) = GenericSite(SVector{length(vec)}(vec))
Base.convert(::Type{GenericSite{N}}, vec::AbstractVector{<:Real}) where N = GenericSite(SVector{N}(vec))
Base.convert(::Type{GenericSite}, tup::Tuple) = GenericSite(SVector(tup))
Base.convert(::Type{GenericSite{N}}, tup::Tuple) where N = GenericSite(SVector{N}(tup))
Base.convert(::Type{GenericSite}, s::AbstractSite) = GenericSite(s.coords)
Base.convert(::Type{GenericSite{N}}, s::AbstractSite{N}) where N = GenericSite(s.coords)

sitekey(s::GenericSite) = round(Int, s.coords[1])

destination(tr::Translation, s::GenericSite) = GenericSite(s.coords + tr.R)

"""
    GenericLattice{SiteT}

A generic lattice of `SiteT` sites.

## Example
```jldoctest
julia> using LatticeModels

julia> l = GenericLattice{2}()
0-site GenericLattice{GenericSite{2}} in 2D space

julia> push!(l, GenericSite(0, 0))  # Add a site at (0, 0)
1-site GenericLattice{GenericSite{2}} in 2D space:
  Site at [0.0, 0.0]

julia> push!(l, (0, 1))             # Add a site at (0, 1)
2-site GenericLattice{GenericSite{2}} in 2D space:
  Site at [0.0, 0.0]
  Site at [0.0, 1.0]

julia> push!(l, [1, 0])             # Add a site at (1, 0)
3-site GenericLattice{GenericSite{2}} in 2D space:
  Site at [0.0, 0.0]
  Site at [0.0, 1.0]
  Site at [1.0, 0.0]

julia> l[2]
2-dim GenericSite{2} at [0.0, 1.0]
```
"""
struct GenericLattice{SiteT} <: AbstractLattice{SiteT}
    sites::Vector{SiteT}
    function GenericLattice(sites::Vector{SiteT}) where SiteT<:AbstractSite{N} where N
        issorted(sites) || throw(ArgumentError("Sites must be sorted"))
        new{SiteT}(sites)
    end
end
GenericLattice(::Vector{<:AbstractSite}) =
    throw(ArgumentError("Sites must be of same dimensionality"))

"""
    GenericLattice{N}()

Constructs an empty `N`-dimensional `GenericLattice` of `GenericSite`s.
"""
GenericLattice{N}() where N = GenericLattice(GenericSite{N}[])
"""
    GenericLattice{SiteType}()

Constructs an empty `GenericLattice` of `SiteType` sites.
"""
GenericLattice{SiteT}() where SiteT<:AbstractSite = GenericLattice(SiteT[])

"""
    GenericLattice(lat)

Constructs a `GenericLattice` from some other lattice `lat`.
"""
GenericLattice(lat::AbstractLattice) = GenericLattice(collect(lat))
GenericLattice(iterable) = GenericLattice([convert(GenericSite, x) for x in iterable])

Base.length(l::GenericLattice) = length(l.sites)
Base.getindex(l::GenericLattice, i::Int) = l.sites[i]
Base.getindex(l::GenericLattice, is::AbstractVector{Int}) = GenericLattice(l.sites[is])
Base.@propagate_inbounds function site_index(l::GenericLattice{SiteT}, site::SiteT, range) where SiteT
    i = searchsortedfirst(@view(l.sites[range]), site) + first(range) - 1
    i in range || return nothing
    return l.sites[i] == site ? i : nothing
end

Base.emptymutable(::GenericLattice, ::Type{SiteT}) where SiteT = GenericLattice(SiteT[])
Base.copymutable(l::GenericLattice) = GenericLattice(copy(l.sites))
function Base.push!(l::GenericLattice{SiteT}, site::SiteT) where SiteT
    i = searchsortedfirst(l.sites, site)
    i ≤ length(l) && l.sites[i] == site && return l
    insert!(l.sites, i, site)
    return l
end
Base.deleteat!(l::GenericLattice, inds) = deleteat!(l.sites, inds)
