using SparseArrays, FillArrays, StaticArrays

abstract type AbstractBonds{LatticeT} end
lattice(bonds::AbstractBonds) = bonds.lat
dims(bonds::AbstractBonds) = dims(lattice(bonds))
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

function Base.iterate(dests::DestinationsIterator, lat_state...)
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
    lat::LT
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

function isdestination(am::AdjacencyMatrix, site1::AbstractSite, site2::AbstractSite)
    i = site_index(am.lat, site1)
    i === nothing && return false
    j = site_index(am.lat, site2)
    j === nothing && return false
    return i < j && am.mat[i, j]
end

nhops(am::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(am.lat, am.mat^n .!= 0)

function adjacency_matrix(bonds::AbstractBonds, more_bonds::AbstractBonds...)
    l = lattice(bonds)
    foreach(more_bonds) do b
        check_samelattice(l, lattice(b))
    end
    Is = Int[]
    Js = Int[]
    for (i, site) in enumerate(l)
        for adj in tuple(bonds, more_bonds...)
            for site2 in destination_sites(adj, site)
                j = site_index(l, site2)
                j === nothing && continue
                push!(Is, i)
                push!(Js, j)
            end
        end
    end
    return AdjacencyMatrix(l, sparse(Is, Js, Fill(true, length(Is)), length(l), length(l), (i,j)->j))
end

abstract type OneToOneBonds{LT} <: AbstractBonds{LT} end
struct UndefinedLattice <: AbstractLattice{NoSite} end
Base.iterate(::UndefinedLattice) = nothing
Base.length(::UndefinedLattice) = 0

function destination(bonds::OneToOneBonds, site::AbstractSite)
    dest = findfirst(site2 -> isdestination(bonds, site, site2), lattice(bonds))
    return dest === nothing ? NoSite() : dest
end
destination_sites(bonds::OneToOneBonds, site::AbstractSite) = (destination(bonds, site),)
function apply_lattice(bonds::AbstractBonds{AbstractLattice}, l::AbstractLattice)
    check_samelattice(l, lattice(bonds))
    return bonds
end
adjacency_matrix(l::AbstractLattice, bonds::AbstractBonds...) =
    adjacency_matrix(Tuple(apply_lattice(b, l) for b in bonds))

Base.:(+)(site::AbstractSite, bonds::OneToOneBonds) = destination(bonds, site)
Base.:(+)(::AbstractSite, ::OneToOneBonds{UndefinedLattice}) =
    throw(ArgumentError("Using a `Bonds`-type object on undefined lattice is allowed only in `build_operator`. Please define the lattice."))
Base.:(-)(bonds::OneToOneBonds) = Base.inv(bonds)
Base.:(-)(site::AbstractSite, bonds::OneToOneBonds) = destination(Base.inv(bonds), site)

struct SpatialShift{LT, N} <: OneToOneBonds{LT}
    lat::LT
    R::SVector{N, Float64}
    function SpatialShift(latt::LT, R::AbstractVector{<:Number}) where
            {N, LT<:AbstractLattice{<:AbstractSite{N}}}
        @check_size R N
        new{LT, N}(latt, SVector{N}(R))
    end
end
function SpatialShift(R::AbstractVector{<:Number})
    n = length(R)
    new{UndefinedLattice, n}(UndefinedLattice(), SVector{n}(R))
end
apply_lattice(sh::SpatialShift{UndefinedLattice}, l::AbstractLattice) =
    SpatialShift(l, sh.R)
dims(sh::SpatialShift{UndefinedLattice, N}) where N = N

isdestination(sh::SpatialShift, site1::AbstractSite, site2::AbstractSite) =
    isapprox(site2.coords - site1.coords, sh.R, atol=√eps())
Base.inv(sh::SpatialShift) = SpatialShift(sh.lat, -sh.R)

Base.summary(io::IO, sh::SpatialShift) = print(io, "Spatial shift with vector R = $(sh.R)")
function Base.show(io::IO, mime::MIME"text/plain", sh::SpatialShift)
    summary(io, sh)
    if !(sh.lat isa UndefinedLattice)
        print(io, "\n on ")
        show(io, mime, sh.lat)
    end
end
