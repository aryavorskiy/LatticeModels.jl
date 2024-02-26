using SparseArrays, FillArrays, StaticArrays, Printf

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
    throw(ArgumentError(
        sprint(show, "text/plain", any, context=(:compact=>true)) *
        " cannot be interpreted as bonds on " *
        sprint(show, "text/plain", l, context=(:compact=>true))))
function adapt_bonds(bonds::AbstractBonds{<:AbstractLattice}, l::AbstractLattice)
    check_samelattice(l, lattice(bonds))
    return bonds
end

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

function Base.show(io::IO, mime::MIME"text/plain", am::AdjacencyMatrix)
    indent = getindent(io)
    print(io, indent, "Adjacency matrix on ")
    summary(io, am.lat)
    requires_compact(io) && return
    print(io, "\n", indent, "Values in a ")
    show(io, mime, am.mat)
end

function adjacentsites(am::AdjacencyMatrix, site::AbstractSite)
    SiteT = eltype(lattice(am))
    rs = resolve_site(am.lat, site)
    rs === nothing && return SiteT[]
    return [rs2.site for rs2 in adjacentsites(am, rs)]
end
function adjacentsites(am::AdjacencyMatrix, rs::ResolvedSite)
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
    DirectedBonds{LT} <: AbstractBonds{LT}

An abstract type for bonds on some lattice that have a direction.

## Methods for subtypes to implement
- `lattice(bonds::DirectionalBonds)`: Returns the lattice where the bonds are defined.
- `destinations(bonds::DirectionalBonds, site::AbstractSite)`: Returns the sites where the
`site` is connected to, accounting for the direction of the bonds.
"""
abstract type DirectedBonds{LT} <: AbstractBonds{LT} end
function destinations end
function isadjacent(bonds::DirectedBonds, site1::AbstractSite, site2::AbstractSite)
    return site2 in destinations(bonds, site1) || site1 in destinations(bonds, site2)
end
function destination(db::DirectedBonds, site::AbstractSite)
    counter = 0
    ret = NoSite()
    for dest in destinations(db, site)
        if dest !== NoSite()
            ret = dest
            counter += 1
            counter > 1 &&
                throw(ArgumentError("The site $site has more than one destination."))
        end
    end
    return ret
end
Base.inv(::DirectedBonds) = throw(ArgumentError("Inverse of the translation is not defined."))
adjacentsites(bonds::DirectedBonds, site::AbstractSite) =
    Base.Iterators.flatten((destinations(bonds, site), destinations(inv(bonds), site)))

function Base.iterate(bonds::DirectedBonds)
    isempty(lattice(bonds)) && return nothing
    l = lattice(bonds)
    isempty(l) && return nothing
    targets = destinations(bonds, first(l))
    jst = iterate(targets)
    return iterate(bonds, (1, targets, jst))
end

@inline function Base.iterate(bonds::DirectedBonds, state)
    i, targets, jst = state
    l = lattice(bonds)
    i > length(l) && return nothing
    if jst === nothing
        i += 1
        i > length(l) && return nothing
        targets = destinations(bonds, l[i])
        jst = iterate(targets)
    end
    jst === nothing && return iterate(bonds, (i, targets, jst))

    rs = ResolvedSite(l[i], i)
    rs2 = resolve_site(l, jst[1])
    jst = iterate(targets, jst[2])

    rs2 === nothing && return iterate(bonds, (i, targets, jst))
    return rs => rs2, (i, targets, jst)
end

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
abstract type AbstractTranslation{LT} <: DirectedBonds{LT} end
destinations(bonds::AbstractTranslation, site::AbstractSite) = (destination(bonds, site),)
@inline function Base.iterate(bonds::AbstractTranslation, i = 1)
    l = lattice(bonds)
    i > length(l) && return nothing
    dest = destination(bonds, l[i])
    s2 = resolve_site(l, dest)
    s2 === nothing && return iterate(bonds, i + 1)
    return ResolvedSite(l[i], i) => s2, i + 1
end
adjacentsites(bonds::AbstractTranslation, site::AbstractSite) =
    (destination(bonds, site), destination(inv(bonds), site))

Base.:(+)(site::AbstractSite, bonds::AbstractTranslation) = destination(bonds, site)
Base.:(+)(::AbstractSite, ::AbstractTranslation{UndefinedLattice}) =
    throw(ArgumentError("Using a `AbstractBonds`-type object on undefined lattice is allowed only in `construct_operator`. Please define the lattice."))
Base.:(-)(bonds::AbstractTranslation) = inv(bonds)
Base.:(-)(site::AbstractSite, bonds::AbstractTranslation) = destination(inv(bonds), site)

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

Base.summary(io::IO, sh::Translation) = print(io, "Translation by ", sh.R)
function Base.show(io::IO, ::MIME"text/plain", sh::Translation)
    indent = getindent(io)
    print(io, indent)
    summary(io, sh)
    requires_compact(io) && return
    if !(sh.lat isa UndefinedLattice)
        print(io, "\n", indent, " on ")
        summary(io, sh.lat)
    end
end

function translate_lattice(l::AbstractLattice, tr::AbstractTranslation)
    e = Base.emptymutable(l, eltype(l))
    ntr = adapt_bonds(tr, l)
    for site in l
        push!(e, destination(ntr, site))
    end
    return e
end
Base.:(+)(l::AbstractLattice, tr::AbstractTranslation) = translate_lattice(l, tr)
