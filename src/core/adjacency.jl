using SparseArrays, FillArrays, StaticArrays

abstract type AbstractBonds{LatticeT} end
lattice(bonds::AbstractBonds) = bonds.lat
dims(bonds::AbstractBonds) = dims(lattice(bonds))
function isadjacent end
isadjacent(bonds::AbstractBonds, s1::ResolvedSite, s2::ResolvedSite) =
    isadjacent(bonds, s1.site, s2.site)
Base.getindex(bonds::AbstractBonds, site1::AbstractSite, site2::AbstractSite) =
    isadjacent(bonds, site1, site2)

function Base.iterate(bonds::AbstractBonds, ind_pair = 1 => 1)
    l = lattice(bonds)
    i, j = ind_pair
    j += 1
    if j > length(l)
        i += 1
        j = i + 1
        j > length(l) && return nothing
    end
    rs1 = ResolvedSite(l[i], i)
    rs2 = ResolvedSite(l[j], j)
    if isadjacent(bonds, rs1, rs2)
        return rs1 => rs2, i => j
    else
        return iterate(bonds, i => j)
    end
end

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
function apply_lattice(b::AdjacencyMatrix, l::AbstractLattice)
    if l == lattice(b)
        return b
    else
        inds = Int[]
        ext_inds = Int[]
        for (i, site) in enumerate(lattice(b))
            j = site_index(l, site)
            if j !== nothing
                push!(inds, i)
                push!(ext_inds, j)
            end
        end
        new_mat = spzeros(Bool, length(l), length(l))
        new_mat[ext_inds, ext_inds] = b.mat[inds, inds]
        return AdjacencyMatrix(l, new_mat)
    end
end

function isadjacent(am::AdjacencyMatrix, site1::AbstractSite, site2::AbstractSite)
    rs1 = resolve_site(am.lat, site1)
    rs1 === nothing && return false
    rs2 = resolve_site(am.lat, site2)
    rs2 === nothing && return false
    return isadjacent(am, rs1, rs2)
end
isadjacent(am::AdjacencyMatrix, s1::ResolvedSite, s2::ResolvedSite) = am.mat[s1.index, s2.index]

nhops(am::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(am.lat, am.mat^n .!= 0)
function Base.union(am::AdjacencyMatrix, ams::AdjacencyMatrix...)
    l = lattice(am)
    for _am in ams
        check_samelattice(l, lattice(_am))
    end
    return AdjacencyMatrix(l, .|(am.mat, getfield.(ams, :mat)...))
end

function adjacency_matrix(bonds::AbstractBonds, more_bonds::AbstractBonds...)
    l = lattice(bonds)
    foreach(more_bonds) do b
        check_samelattice(l, lattice(b))
    end
    Is = Int[]
    Js = Int[]
    for adj in tuple(bonds, more_bonds...)
        for (s1, s2) in adj
            push!(Is, s1.index)
            push!(Js, s2.index)
        end
    end
    return AdjacencyMatrix(l, sparse(Is, Js, Fill(true, length(Is)), length(l), length(l), (i,j)->j))
end
adjacency_matrix(l::AbstractLattice, bonds::AbstractBonds...) =
    adjacency_matrix((apply_lattice(b, l) for b in bonds)...)

struct UndefinedLattice <: AbstractLattice{NoSite} end
Base.iterate(::UndefinedLattice) = nothing
Base.length(::UndefinedLattice) = 0

abstract type AbstractTranslation{LT} <: AbstractBonds{LT} end
isadjacent(bonds::AbstractTranslation, site1::AbstractSite, site2::AbstractSite) =
    site2 === destination(bonds, site1) || site1 === destination(bonds, site2)
@inline function Base.iterate(bonds::AbstractTranslation, i = 1)
    l = lattice(bonds)
    i > length(l) && return nothing
    dest = destination(bonds, l[i])
    s2 = resolve_site(l, dest)
    s2 === nothing && return iterate(bonds, i + 1)
    return ResolvedSite(l[i], i) => s2, i + 1
end
function apply_lattice(bonds::AbstractBonds{<:AbstractLattice}, l::AbstractLattice)
    check_samelattice(l, lattice(bonds))
    return bonds
end

Base.:(+)(site::AbstractSite, bonds::AbstractTranslation) = destination(bonds, site)
Base.:(+)(::AbstractSite, ::AbstractTranslation{UndefinedLattice}) =
    throw(ArgumentError("Using a `AbstractBonds`-type object on undefined lattice is allowed only in `build_operator`. Please define the lattice."))
Base.inv(::AbstractTranslation) = throw(ArgumentError("Inverse of the translation is not defined."))
Base.:(-)(bonds::AbstractTranslation) = Base.inv(bonds)
Base.:(-)(site::AbstractSite, bonds::AbstractTranslation) = destination(Base.inv(bonds), site)

struct SpatialShift{LT, N} <: AbstractTranslation{LT}
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
apply_lattice(bonds::SpatialShift, l::AbstractLattice) =
    SpatialShift(l, bonds.R)
function destination(sh::SpatialShift, site::AbstractSite)
    for dest in lattice(sh)
        if isapprox(site.coords + sh.R, dest.coords, atol=√eps())
            return dest
        end
    end
    return NoSite()
end
dims(::SpatialShift{UndefinedLattice, N}) where N = N

isadjacent(sh::SpatialShift, site1::AbstractSite, site2::AbstractSite) =
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
