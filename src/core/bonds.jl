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

abstract type BondsIterateStyle end
struct IterateDirect <: BondsIterateStyle end
struct IterateNearest <: BondsIterateStyle end
iterstyle(::AbstractBonds) = IterateDirect()
@inline function _destination_sites(::IterateDirect, bonds::AbstractBonds, site::AbstractSite)
    l = lattice(bonds)
    return (site2 for site2 in l
        if isadjacent(bonds, site, site2) && site2 > site)
end
@inline function _destination_sites(::IterateNearest, bonds::AbstractBonds, site::AbstractSite)
    l = lattice(bonds)
    return (translate_to_nearest(l, site, site2) for site2 in l
        if isadjacent(bonds, site, translate_to_nearest(l, site, site2))
            && site2 > site)
end
@inline _destination_sites(bonds::AbstractBonds, site::AbstractSite) =
    _destination_sites(iterstyle(bonds), bonds, site)

abstract type BondsResolveStyle end
struct ResolveIndices <: BondsResolveStyle end
struct ResolveSites <: BondsResolveStyle end
resolvestyle(::AbstractBonds) = ResolveSites()
function _resolve_all(::IterateDirect, bonds::AbstractBonds, site, dests)
    l = lattice(bonds)
    return (resolve_site(l, dest) for dest in dests)
end
function _resolve_all(::IterateNearest, bonds::AbstractBonds, site, dests)
    l = lattice(bonds)
    return (resolve_site(l, translate_to_nearest(l, site, dest)) for dest in dests)
end
_destinations(::ResolveIndices, bonds::AbstractBonds, rs::ResolvedSite) =
    _resolve_all(iterstyle(bonds), bonds, rs.site, _destination_inds(bonds, rs.index))
_destinations(::ResolveSites, bonds::AbstractBonds, rs::ResolvedSite) =
    _resolve_all(iterstyle(bonds), bonds, rs.site, _destination_sites(bonds, rs.site))
_destinations(bonds::AbstractBonds, rs::ResolvedSite) =
    _destinations(resolvestyle(bonds), bonds, rs)

function Base.iterate(bonds::AbstractBonds)
    l = lattice(bonds)
    rs = resolve_site(l, 1)
    rs === nothing && return nothing
    targets = _destinations(bonds, rs)
    jst = iterate(targets)
    return iterate(bonds, (1, targets, jst))
end

@inline function Base.iterate(bonds::AbstractBonds, state)
    i, targets, jst = state
    l = lattice(bonds)
    i > length(l) && return nothing
    if jst === nothing
        i += 1
        i > length(l) && return nothing
        targets = _destinations(bonds, resolve_site(l, i))
        jst = iterate(targets)
    end
    jst === nothing && return iterate(bonds, (i, targets, jst))

    rs = resolve_site(l, i)
    rs2 = jst[1]
    jst = iterate(targets, jst[2])

    rs2 === nothing && return iterate(bonds, (i, targets, jst))
    return rs => rs2, (i, targets, jst)
end

"""
    translate_to_nearest(lat, site1, site2)

Translate `site2` to its equivalent nearest to `site1` in the lattice `lat`, taking the
boundary conditions into account.
"""
translate_to_nearest(::AbstractLattice, ::AbstractSite, site2::AbstractSite) = site2
translate_to_nearest(l::AbstractLattice, site1::AbstractSite, i2::Int) =
    translate_to_nearest(l, site1, l[i2])

"""
    sitedistance([lat, ]site1, site2)
Returns the distance between two sites on the `lat` lattice, taking boundary conditions into account.

# Arguments
- `lat`: The lattice where the sites are defined.
- `site1` and `site2`: The sites to measure the distance between.
"""
sitedistance(::AbstractLattice, site1::AbstractSite, site2::AbstractSite) =
    norm(site1.coords - site2.coords)
sitedistance(site1, site2) = sitedistance(UndefinedLattice(), site1, site2)
sitedistance(l::LatticeWithMetadata, site1::AbstractSite, site2::AbstractSite) =
    norm(translate_to_nearest(l, site1, site2).coords - site1.coords)

"""
    SiteDistance(f, lat)

A bonds type that connects sites based on the distance between them.

# Arguments
- `f`: A function that takes a distance and returns if the distance is allowed.
- `lat`: The lattice where the bonds are defined.
"""
struct SiteDistance{LT, FT} <: AbstractBonds{LT}
    f::FT
    lat::LT
end
SiteDistance(f) = SiteDistance(f, UndefinedLattice())
iterstyle(::SiteDistance) = IterateNearest()
isadjacent(bonds::SiteDistance{<:Any, <:Function}, s1::AbstractSite, s2::AbstractSite) =
    bonds.f(sitedistance(bonds.lat, s1, s2))
isadjacent(bonds::SiteDistance, s1::AbstractSite, s2::AbstractSite) =
    sitedistance(bonds.lat, s1, s2) in bonds.f
adapt_bonds(bonds::SiteDistance, lat::AbstractLattice) = SiteDistance(bonds.f, lat)
function Base.show(io::IO, ::MIME"text/plain", sd::SiteDistance)
    indent = getindent(io)
    print(io, indent, "SiteDistance(")
    show(io, sd.f)
    if sd.lat !== UndefinedLattice()
        print(io, ", ")
        summary(io, sd.lat)
    end
    print(io, ")")
end

"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

---
    AdjacencyMatrix(lat[, mat])

Construct an adjacency matrix from the `mat` matrix on the `lat` lattice.

If `mat` is not provided, it is assumed to be a zero matrix.

## Example
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(2, 2);

julia> a = AdjacencyMatrix(l)
Adjacency matrix on 4-site SquareLattice in 2D space
Values in a 4×4 SparseArrays.SparseMatrixCSC{Bool, Int64} with 0 stored entries:
 ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅

julia> site1, site2, site3, site4 = l;

julia> a[site1, site2] = a[site2, site4] = a[site3, site4] = true;

julia> a
Adjacency matrix on 4-site SquareLattice in 2D space
Values in a 4×4 SparseArrays.SparseMatrixCSC{Bool, Int64} with 6 stored entries:
 ⋅  1  ⋅  ⋅
 1  ⋅  ⋅  1
 ⋅  ⋅  ⋅  1
 ⋅  1  1  ⋅
```
"""
struct AdjacencyMatrix{LT,MT,ST} <: AbstractBonds{LT}
    lat::LT
    mat::MT
    iterstyle::ST
    function AdjacencyMatrix(lat::LT, mat::MT, style::ST=IterateDirect()) where {LT<:AbstractLattice,MT<:AbstractMatrix{Bool}, ST<:BondsIterateStyle}
        @check_size mat :square
        @check_size lat size(mat, 1)
        eye = spdiagm(Fill(true, length(lat)))
        new{LT,MT,ST}(lat, dropzeros((mat .| transpose(mat)) .& .!eye), style)
    end
end
function AdjacencyMatrix(l::AbstractLattice, style::ST=IterateDirect()) where {ST<:BondsIterateStyle}
    AdjacencyMatrix(l, spzeros(Bool, length(l), length(l)), style)
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
        return AdjacencyMatrix(l, new_mat, b.iterstyle)
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

iterstyle(adj::AdjacencyMatrix) = adj.iterstyle
resolvestyle(::AdjacencyMatrix) = ResolveIndices()
@inline _destination_inds(am::AdjacencyMatrix, i::Int) =
    (j for j in i + 1:length(lattice(am)) if am.mat[i, j])
@inline function _destination_inds(am::AdjacencyMatrix{<:Any,<:SparseMatrixCSC}, i::Int)
    p1 = am.mat.colptr[i]
    p2 = am.mat.colptr[i + 1]
    st = findfirst(>(i), @view(am.mat.rowval[p1:p2 - 1]))
    st === nothing && return ()
    v = @view(am.mat.rowval[p1 + st - 1:p2 - 1])
    return (v[k] for k in eachindex(v) if am.mat.nzval[p1 + st - 1 + k - 1])
end

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

"""
    adjacentsites(bonds, site)

Returns the sites that are connected to `site` by the `bonds`.
"""
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

"""
    AdjacencyMatrix([lat, ]bonds...)

Constructs an adjacency matrix from the `bonds`. If `lat` is not provided, it is inferred
from the `bonds`.
"""
function AdjacencyMatrix(bonds::AbstractBonds, more_bonds::AbstractBonds...)
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
    mat = sparse(Is, Js, Fill(true, length(Is)), length(l), length(l), (i,j)->j)
    return AdjacencyMatrix(l, mat)
end
AdjacencyMatrix(l::AbstractLattice, bonds::AbstractBonds...) =
AdjacencyMatrix((adapt_bonds(b, l) for b in bonds)...)

"""
    AdjacencyMatrix(f, lat)

Constructs an adjacency matrix from the function `f` that returns if the sites are connected
on the `lat` lattice.
"""
function AdjacencyMatrix(f::Function, l::AbstractLattice)
    Is = Int[]
    Js = Int[]
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            j > i || continue
            if f(site1, site2)
                push!(Is, i, j)
                push!(Js, j, i)
            end
        end
    end
    mat = sparse(Is, Js, Fill(true, length(Is)), length(l), length(l))
    return AdjacencyMatrix(l, mat)
end

"""
    UndefinedLattice

A lattice that is not defined.
The bonds can be 'defined' on it in context where the lattice is already defined before,
e. g. in `construct_operator`.
"""
struct UndefinedLattice <: AbstractLattice{NoSite} end
Base.iterate(::UndefinedLattice) = nothing
Base.length(::UndefinedLattice) = 0
Base.copymutable(::UndefinedLattice) = UndefinedLattice()
Base.emptymutable(::UndefinedLattice) = UndefinedLattice()

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
_destination_sites(bonds::DirectedBonds, site::AbstractSite) =
    destinations(bonds, site)
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
destination(tr::AbstractTranslation, site::AbstractSite) = throw(MethodError(destination, (tr, site)))
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

Base.:(+)(site::AbstractSite, bonds::DirectedBonds) = destination(bonds, site)
Base.:(-)(bonds::DirectedBonds) = inv(bonds)
Base.:(-)(site::AbstractSite, bonds::DirectedBonds) = destination(inv(bonds), site)

"""
    Translation <: AbstractTranslation

A spatial translation on some lattice.

## Fields
- `lat`: The lattice where the translations are defined.
- `R`: The vector of the translation.

## Example
```jldoctest
julia> using LatticeModels

julia> gl = GenericLattice([(1, 1), (1, 2), (2, 1), (2, 2)])
4-site GenericLattice{GenericSite{2}} in 2D space:
  Site at [1.0, 1.0]
  Site at [1.0, 2.0]
  Site at [2.0, 1.0]
  Site at [2.0, 2.0]

julia> tr = Translation(gl, [1, 0])     # Translation by [1, 0]
Translation by [1.0, 0.0]
 on 4-site GenericLattice{GenericSite{2}} in 2D space

julia> site1 = gl[!, x = 1, y = 1]      # Site at [1, 1]
2-dim GenericSite{2} at [1.0, 1.0]

julia> site1 + tr                       # Translated site
2-dim GenericSite{2} at [2.0, 1.0]

julia> site1 - tr                       # Inverse translation
2-dim GenericSite{2} at [0.0, 1.0]
```
"""
struct Translation{LT, N} <: AbstractTranslation{LT}
    lat::LT
    R::SVector{N, Float64}
    function Translation(lat::LT, R::AbstractVector{<:Number}) where
            {N, LT<:AbstractLattice{<:AbstractSite{N}}}
        @check_size R N
        new{LT, N}(lat, SVector{N}(R))
    end
    function Translation(R::AbstractVector{<:Number})
        n = length(R)
        new{UndefinedLattice, n}(UndefinedLattice(), SVector{n}(R))
    end
end
Translation(::UndefinedLattice, R::AbstractVector{<:Number}) = Translation(R)
adapt_bonds(bonds::Translation, l::AbstractLattice) = Translation(l, bonds.R)
adapt_bonds(bonds::Translation, ::UndefinedLattice) = Translation(bonds.R)
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

function translate_lattice(l::AbstractLattice, tr::DirectedBonds)
    e = Base.emptymutable(l, eltype(l))
    ntr = adapt_bonds(tr, l)
    for site in l
        push!(e, destination(ntr, site))
    end
    return e
end
Base.:(+)(l::AbstractLattice, tr::DirectedBonds) = translate_lattice(l, tr)
