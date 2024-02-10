using SparseArrays, FillArrays, StaticArrays

abstract type AbstractBonds{LatticeT} end
function lattice end
function isdestination end
Base.getindex(bonds::AbstractBonds, site1::AbstractSite, site2::AbstractSite) =
    isdestination(bonds, site1, site2)

struct DestinationsIterator{AT, ST}
    bonds::AT
    site::ST
    function DestinationsIterator(bonds::AT, site::ST) where
            {ST<:AbstractSite, AT<:AbstractBonds{<:AbstractLattice{ST}}}
        new{AT, ST}(bonds, site)
    end
end

function iterate(dests::DestinationsIterator, lat_state...)
    p = iterate(lattice(dests.bonds), lat_state...)
    p === nothing && return nothing
    site, new_state = p
    if isdestination(dests.bonds, dests.site, site) && site != dests.site
        return site, new_state
    else
        return iterate(dests, new_state)
    end
end
destination_sites(bonds::AbstractBonds, site::AbstractSite) = DestinationsIterator(bonds, site)

"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`AdjacencyMatrix`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `AdjacencyMatrix` which connects sites that were connected by `≤n` bonds of the previous `AdjacencyMatrix`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct AdjacencyMatrix{LT,MT} <: AbstractBonds{LT}
    lattice::LT
    mat::MT
    function AdjacencyMatrix(l::LT, connectivity_matrix::MT) where {LT<:AbstractLattice,MT<:AbstractMatrix{Bool}}
        @check_size connectivity_matrix :square
        @check_size l size(connectivity_matrix, 1)
        new{LT,MT}(l,
            (connectivity_matrix .| transpose(connectivity_matrix)) .& .!Eye{Bool}(length(l)))
    end
    function AdjacencyMatrix(l::AbstractLattice)
        AdjacencyMatrix(l, spzeros(length(l), length(l)))
    end
end

lattice(am::AdjacencyMatrix) = am.lattice

function isdestination(am::AdjacencyMatrix, site1::AbstractSite, site2::AbstractSite)
    i = site_index(am.lattice, site1)
    i === nothing && return false
    j = site_index(am.lattice, site2)
    j === nothing && return false
    return i < j && am.bmat[i, j]
end

nhops(am::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(am.sites, am.bmat^n .!= 0)

function AdjacencyMatrix(adj::AbstractBonds)
    l = lattice(adj)
    Is = Int[]
    Js = Int[]
    for (i, site) in enumerate(l)
        for site2 in destination_sites(adj, site)
            j = site_index(l, site2)
            j === nothing && continue
            push!(Is, i)
            push!(Js, j)
        end
    end
    return AdjacencyMatrix(l, sparse(Is, Js, Fill(1, length(Is))))
end

abstract type OneToOneBonds{LT} <: AbstractBonds{LT} end
function destination(bonds::OneToOneBonds, site::AbstractSite)
    dest = findfirst(site2 -> isdestination(bonds, site, site2), lattice(bonds))
    return dest === nothing ? NoSite() : dest
end
destination_sites(bonds::OneToOneBonds, site::AbstractSite) = (destination(bonds, site),)

Base.:(+)(site::AbstractSite, bonds::OneToOneBonds) = destination(bonds, site)
Base.:(-)(bonds::OneToOneBonds) = Base.inv(bonds)
Base.:(-)(site::AbstractSite, bonds::OneToOneBonds) = destination(Base.inv(bonds), site)

struct SpatialShift{LT, N}
    latt::LT
    translation_vec::SVector{Float64, N}
    function SpatialShift(latt::LT, tr::SVector{<:Number,N}) where
            {N, LT<:AbstractLattice{<:AbstractSite{N}}}
        new{LT, N}(latt, tr)
    end
end

isdestination(sh::SpatialShift, site1::AbstractSite, site2::AbstractSite) =
    isapprox(site2.coords - site1.coords, sh.translation_vec, atol=√eps())
Base.inv(sh::SpatialShift) = SpatialShift(sh.latt, -sh.translation_vec)
