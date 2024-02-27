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
Base.:(==)(s1::GenericSite, s2::GenericSite) = s1.coords == s2.coords
Base.isless(s1::GenericSite, s2::GenericSite) = isless(s1.coords, s2.coords)

Base.convert(::Type{GenericSite}, vec::AbstractVector{<:Real}) = GenericSite(SVector{length(vec)}(vec))
Base.convert(::Type{GenericSite{N}}, vec::AbstractVector{<:Real}) where N = GenericSite(SVector{N}(vec))
Base.convert(::Type{GenericSite}, tup::Tuple) = GenericSite(SVector(tup))
Base.convert(::Type{GenericSite{N}}, tup::Tuple) where N = GenericSite(SVector{N}(tup))

sitekey(s::GenericSite) = round(Int, s.coords[1])

"""
    GenericLattice{SiteT}

A generic lattice of `SiteT` sites.
"""
struct GenericLattice{SiteT} <: AbstractLattice{SiteT}
    sites::Vector{SiteT}
end
GenericLattice{N}() where N = GenericLattice(GenericSite{N}[])
GenericLattice(lat::AbstractLattice) = GenericLattice(collect(lat))

Base.length(l::GenericLattice) = length(l.sites)
Base.getindex(l::GenericLattice, i::Int) = l.sites[i]
Base.getindex(l::GenericLattice, is::AbstractVector{Int}) = GenericLattice(l.sites[is])
Base.@propagate_inbounds function site_index(l::GenericLattice{SiteT}, site::SiteT, range) where SiteT
    i = searchsortedfirst(@view(l.sites[range]), site) + first(range) - 1
    i in range || return nothing
    return l.sites[i] == site ? i : nothing
end
site_index(l::GenericLattice{SiteT}, site::SiteT) where SiteT<:AbstractSite = findfirst(==(site), l.sites)

Base.emptymutable(::GenericLattice, ::Type{SiteT}) where SiteT = GenericLattice(SiteT[])
Base.copymutable(l::GenericLattice) = GenericLattice(copy(l.sites))
function Base.push!(l::GenericLattice{SiteT}, site::SiteT) where SiteT
    i = searchsortedfirst(l.sites, site)
    i ≤ length(l) && l.sites[i] == site && return l
    insert!(l.sites, i, site)
    return l
end
Base.deleteat!(l::GenericLattice, inds) = deleteat!(l.sites, inds)
