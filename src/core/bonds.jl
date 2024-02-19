using SparseArrays, FillArrays, StaticArrays

const SingleBond{LT<:AbstractSite} = Pair{LT, LT}

"""
    AbstractBonds{LT}

An abstract type for bonds on some lattice.

## Methods for subtypes to implement
- `lattice(bonds::AbstractBonds)`: Returns the lattice where the bonds are defined.
- `isadjacent(bonds::AbstractBonds, site1::AbstractSite, site2::AbstractSite)`:
    Returns if the sites are connected by the bonds.

## Optional methods for subtypes to implement
- `adapt_bonds(bonds::AbstractBonds, l::AbstractLattice)`
"""
abstract type AbstractBonds{LatticeT} end
lattice(bonds::AbstractBonds) = bonds.lat
dims(bonds::AbstractBonds) = dims(lattice(bonds))
function isadjacent end
isadjacent(bonds::AbstractBonds, s1::ResolvedSite, s2::ResolvedSite) =
    isadjacent(bonds, s1.site, s2.site)
Base.getindex(bonds::AbstractBonds, site1::AbstractSite, site2::AbstractSite) =
    isadjacent(bonds, site1, site2)

"""
    adapt_bonds(bonds, lat)

Adapt the bonds to the lattice `lat`. The output can be a different type of
bonds, more fitting for the concrete type of lattice.
"""
adapt_bonds(any, l::AbstractLattice) =
    throw(ArgumentError(sprint(show, "text/plain", any) * " cannot be interpreted as bonds on " * sprint(show, "text/plain", l)))

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
    site_distance([lat, ]site1, site2)
Returns the distance between two sites on the `l` lattice, taking boundary conditions into account.

# Arguments
- `l::Lattice`: The lattice where the sites are defined.
- `site1::LatticeSite`: The first site.
"""
site_distance(::AbstractLattice, site1::AbstractSite, site2::AbstractSite) =
    norm(site1.coords - site2.coords)
site_distance(site1, site2) = site_distance(UndefinedLattice(), site1, site2)

"""
    SiteDistance{FT}

A bonds type that connects sites based on the distance between them.

# Arguments
- `f::Function`: A function that takes a distance and returns if the distance is allowed.
"""
struct SiteDistance{LT, FT} <: AbstractBonds{LT}
    lat::LT
    f::Function
end

isadjacent(bonds::SiteDistance, s1::AbstractSite, s2::AbstractSite) =
    bonds.f(site_distance(bonds.lat, s1, s2))
adapt_bonds(bonds::SiteDistance, ::AbstractLattice) = bonds

"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.
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
function adapt_bonds(b::AdjacencyMatrix, l::AbstractLattice)
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

function Base.setindex!(b::AdjacencyMatrix, v, site1::AbstractSite, site2::AbstractSite)
    rs1 = resolve_site(b.lat, site1)
    rs1 === nothing && return nothing
    rs2 = resolve_site(b.lat, site2)
    rs2 === nothing && return nothing
    setindex!(b, v, rs1, rs2)
end
function Base.setindex!(b::AdjacencyMatrix, v, s1::ResolvedSite, s2::ResolvedSite)
    b.mat[s1.index, s2.index] = v
    b.mat[s2.index, s1.index] = v
end

function Base.union(am::AdjacencyMatrix, ams::AdjacencyMatrix...)
    l = lattice(am)
    for _am in ams
        check_samelattice(l, lattice(_am))
    end
    return AdjacencyMatrix(l, .|(am.mat, getfield.(ams, :mat)...))
end

function target_sites(am::AdjacencyMatrix, site::AbstractSite)
    SiteT = eltype(lattice(am))
    rs = resolve_site(am.lat, site)
    rs === nothing && return SiteT[]
    return [rs2.site for rs2 in target_sites(am, rs)]
end
function target_sites(am::AdjacencyMatrix, rs::ResolvedSite)
    l = lattice(am)
    return [ResolvedSite(l[i], i) for i in findall(i->am.mat[rs.index, i], eachindex(l))]
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
    adjacency_matrix((adapt_bonds(b, l) for b in bonds)...)

"""
    UndefinedLattice

A lattice that is not defined.
The bonds can be 'defined' on it in context where the lattice is already defined before,
e. g. in `construct_operator`.
"""
struct UndefinedLattice <: AbstractLattice{NoSite} end
Base.iterate(::UndefinedLattice) = nothing
Base.length(::UndefinedLattice) = 0

"""
    AbstractTranslation{LT}

An abstract type for translations on some lattice.

## Methods for subtypes to implement
- `lattice(bonds::AbstractTranslation)`: Returns the lattice where the translations are defined.
- `destination(bonds::AbstractTranslation, site::AbstractSite)`: Returns the site where the `site` is translated to.

## Optional methods for subtypes to implement
- `adapt_bonds(bonds::AbstractTranslation, l::AbstractLattice)`:
    Adapt the translation to the lattice `l`. The output can be a different type of
    translation, more fitting for the concrete type of lattice.
- `inv(bonds::AbstractTranslation)`: Returns the inverse of the translation, if any.
"""
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
function adapt_bonds(bonds::AbstractBonds{<:AbstractLattice}, l::AbstractLattice)
    check_samelattice(l, lattice(bonds))
    return bonds
end

Base.:(+)(site::AbstractSite, bonds::AbstractTranslation) = destination(bonds, site)
Base.:(+)(::AbstractSite, ::AbstractTranslation{UndefinedLattice}) =
    throw(ArgumentError("Using a `AbstractBonds`-type object on undefined lattice is allowed only in `construct_operator`. Please define the lattice."))
Base.inv(::AbstractTranslation) = throw(ArgumentError("Inverse of the translation is not defined."))
Base.:(-)(bonds::AbstractTranslation) = Base.inv(bonds)
Base.:(-)(site::AbstractSite, bonds::AbstractTranslation) = destination(Base.inv(bonds), site)

"""
    Translation <: AbstractTranslation

A spatial translation on some lattice.

## Fields
- `lat`: The lattice where the translations are defined.
- `R`: The vector of the translation.
"""
struct Translation{LT, N} <: AbstractTranslation{LT}
    lat::LT
    R::SVector{N, Float64}
    function Translation(latt::LT, R::AbstractVector{<:Number}) where
            {N, LT<:AbstractLattice{<:AbstractSite{N}}}
        @check_size R N
        new{LT, N}(latt, SVector{N}(R))
    end
    function Translation(R::AbstractVector{<:Number})
        n = length(R)
        new{UndefinedLattice, n}(UndefinedLattice(), SVector{n}(R))
    end
end
adapt_bonds(bonds::Translation, l::AbstractLattice) = Translation(l, bonds.R)
adapt_bonds(bonds::Translation, ::UndefinedLattice) = Translation(bonds.R)
function destination(sh::Translation, site::AbstractSite)
    for dest in lattice(sh)
        if isapprox(site.coords + sh.R, dest.coords, atol=√eps())
            return dest
        end
    end
    return NoSite()
end
dims(::Translation{UndefinedLattice, N}) where N = N

isadjacent(sh::Translation, site1::AbstractSite, site2::AbstractSite) =
    isapprox(site2.coords - site1.coords, sh.R, atol=√eps())
Base.inv(sh::Translation) = Translation(sh.lat, -sh.R)

Base.summary(io::IO, sh::Translation) = print(io, "Spatial shift with vector R = $(sh.R)")
function Base.show(io::IO, mime::MIME"text/plain", sh::Translation)
    summary(io, sh)
    if !(sh.lat isa UndefinedLattice)
        print(io, "\n on ")
        summary(io, sh.lat)
    end
end

"""
    NearestNeighbor{N}

A bonds type that connects sites that are nearest neighbors of order `N` on some lattice.
"""
struct NearestNeighbor{N} <: AbstractBonds{UndefinedLattice}
    function NearestNeighbor(::Val{N}) where {N}
        new{N}()
    end
end
NearestNeighbor(N::Int) = NearestNeighbor(Val(N))

struct DefaultNNBonds{M, TupleT}
    dists::NTuple{M, Float64}
    nnbonds::TupleT
    function DefaultNNBonds(dists::NTuple{M, Float64}, nnbonds::Tuple) where {M}
        length(dists) == length(nnbonds) ||
            throw(ArgumentError("The number of distances and the number of hops must be the same."))
        nnbonds_nolat = Tuple(adapt_bonds(b, UndefinedLattice()) for b in nnbonds)
        new{M, typeof(nnbonds_nolat)}(dists, nnbonds_nolat)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", dnn::DefaultNNBonds)
    io = IOContext(io, :inline => true)
    for i in 1:length(dnn.dists)
        println(io, "\n", dnn.dists[i], " =>")
        show(io, mime, dnn.nnbonds[i])
    end
end

Base.getindex(dnn::DefaultNNBonds, i::Int) = dnn.nnbonds[i]
Base.length(dnn::DefaultNNBonds) = length(dnn.nnbonds)

getnnbonds(l::AbstractLattice) = getparam(l, :nnbonds, DefaultNNBonds((), ()))
setnnbonds(l::AbstractLattice, dnn::DefaultNNBonds) = setparam(l, :nnbonds, dnn)

function adapt_bonds(b::NearestNeighbor{N}, l::LatticeWithParams) where {N}
    default_nnhops = getnnbonds(l)
    if default_nnhops === nothing || N > length(default_nnhops)
        return adapt_bonds(adapt_bonds(b, l.lat), l)
    else
        return adapt_bonds(default_nnhops[N], l)
    end
end

function transform_lattice(l::AbstractLattice, tr::AbstractTranslation)
    e = Base.emptymutable(l, eltype(l))
    ntr = adapt_bonds(tr, l)
    for site in l
        push!(e, destination(ntr, site))
    end
    return e
end
Base.:(+)(l::AbstractLattice, tr::AbstractTranslation) = transform_lattice(l, tr)
