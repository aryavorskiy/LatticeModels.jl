"""
    GenericSite{N}

A generic site in an N-dimensional lattice.
"""
struct GenericSite{N} <: AbstractSite{N}
    coord::SVector{N,Int}
end
Base.:(==)(s1::GenericSite, s2::GenericSite) = s1.coord == s2.coord

"""
    GenericLattice{SiteT}

A generic lattice of `SiteT` sites.
"""
struct GenericLattice{SiteT} <: AbstractLattice{SiteT}
    sites::Vector{SiteT}
end
GenericLattice(lat::AbstractLattice) = GenericLattice(collect(lat))

Base.length(l::GenericLattice) = length(l.sites)
Base.getindex(l::GenericLattice, i::Int) = l.sites[i]
Base.getindex(l::GenericLattice, is::AbstractVector{Int}) = GenericLattice(l.sites[is])
site_index(l::GenericLattice, site::SiteT) = findfirst(==(site), l.sites)

Base.emptymutable(::GenericLattice, ::Type{SiteT}) where SiteT = GenericLattice(SiteT[])
Base.copymutable(l::GenericLattice) = GenericLattice(copy(l.sites))
Base.push!(l::GenericLattice, site::SiteT) = push!(l.sites, site)
Base.deleteat!(l::GenericLattice, inds) = deleteat!(l.sites, inds)
