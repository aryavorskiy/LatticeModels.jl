import Base: getindex, +, -, *, /
using LinearAlgebra

"""
    AbstractCurrents

Supertype for all type representing currents-like values on a lattice.
Subtypes must implement `Base.getindex(::Int, ::Int)` and `lattice` functions.
"""
abstract type AbstractCurrents end

"""
Gets the lattice where the given `AbstractCurrents` object is defined.
"""
lattice(curr::AbstractCurrents) = error("lattice(::$(typeof(curr))) must be explicitly implemented")

Base.getindex(curr::AbstractCurrents, s1::AbstractSite, s2::AbstractSite) =
    curr[site_index(lattice(curr), s1), site_index(lattice(curr), s2)]
Base.getindex(curr::AbstractCurrents, s1::ResolvedSite, s2::ResolvedSite) =
    curr[s1.index, s2.index]

function Base.show(io::IO, mime::MIME"text/plain", curr::AbstractCurrents)
    summary(io, curr)
    print(io, " on ")
    io = IOContext(io, :compact => true)
    show(io, mime, lattice(curr))
end

@inline iszerocurrent(::AbstractCurrents, i::Int, j::Int) = false
@inline iszerocurrent(curr::AbstractCurrents, s1::AbstractSite, s2::AbstractSite) =
    iszerocurrent(curr, site_index(lattice(curr), s1), site_index(lattice(curr), s2))

function Base.length(curr::AbstractCurrents)
    le = length(lattice(curr))
    return le * (le - 1) ÷ 2
end
@inline function Base.iterate(curr::AbstractCurrents)
    if length(lattice(curr)) < 2
        return nothing
    else
        lat = lattice(curr)
        return (lat[1] => lat[2], curr[1, 2]), 1 => 2
    end
end
@inline function Base.iterate(curr::AbstractCurrents, pair::Pair{Int, Int})
    lat = lattice(curr)
    i, j = pair
    j += 1
    if j > length(lat)
        i += 1
        j = i + 1
        j > length(lat) && return nothing
    end
    iszerocurrent(curr, i, j) && return iterate(curr, i => j)
    return (lat[i] => lat[j], curr[i, j]), i => j
end

"""
    SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents

A lazy wrapper for a `Currents` object representing the same currents but on a smaller lattice.
"""
struct SubCurrents{CT} <: AbstractCurrents
    parent_currents::CT
    lattice::AbstractLattice
    indices::Vector{Int}
    function SubCurrents(parent_currents::CT, indices::Vector{Int}) where {CT<:AbstractCurrents}
        lat = lattice(parent_currents)
        new{CT}(parent_currents, lat[indices], indices)
    end
end

Base.getindex(scurr::SubCurrents, i::Int, j::Int) = scurr.parent_currents[scurr.indices[i], scurr.indices[j]]
lattice(scurr::SubCurrents) = scurr.lattice
function Base.getindex(curr::AbstractCurrents, lvm::LatticeValue{Bool})
    check_samesites(curr, lvm)
    indices = findall(lvm.values)
    SubCurrents(curr, indices)
end

function Base.summary(io::IO, scurr::SubCurrents)
    print(io, "SubCurrents of ")
    summary(io, scurr.parent_currents)
end

function _site_indices(l::AbstractLattice, l2::AbstractLattice)
    check_issublattice(l2, l)
    return [site_index(l, site) for site in l2]
end
_site_indices(l::AbstractLattice, site::AbstractSite) = (site_index(l, site),)
function currents_from(curr::AbstractCurrents, src)
    lat = lattice(curr)
    is = _site_indices(lat, src)
    LatticeValue(lat, Float64[j in is ? 0 : sum(curr[i, j] for i in is) for j in eachindex(lat)])
end
function currents_from_to(curr::AbstractCurrents, src, dst=nothing)
    lat = lattice(curr)
    is = _site_indices(lat, src)
    js = dst === nothing ? setdiff(eachindex(lat), is) : _site_indices(lat, dst)
    sum(curr[i, j] for i in is, j in js)
end

"""
    Currents <: AbstractCurrents

A `AbstractCurrents` instance that stores values for all currents explicitly.
"""
struct Currents{MT,LT} <: AbstractCurrents
    lat::LT
    currents::MT
    function Currents(l::LT, curs::MT) where {MT<:AbstractMatrix, LT}
        @check_size curs :square
        @check_size l size(curs, 1)
        new{MT, LT}(l, curs)
    end
end
Currents(l::AbstractLattice) = Currents(l, spzeros(length(l), length(l)))

Base.convert(::Type{Currents}, curr::AbstractCurrents) = Currents(curr)
Base.copy(mc::Currents) = Currents(lattice(mc), copy(mc.currents))
Base.zero(mc::Currents) = Currents(lattice(mc), zero(mc.currents))
lattice(mcurr::Currents) = mcurr.lat

Base.:(==)(c1::Currents, c2::Currents) =
    c1.lat == c2.lat && c1.currents == c2.currents
Base.isapprox(c1::Currents, c2::Currents; kw...) =
    c1.lat == c2.lat && isapprox(c1.currents, c2.currents; kw...)

Base.getindex(mcurr::Currents, i::Int, j::Int) = mcurr.currents[i, j]
function Base.setindex!(curr::Currents, rhs, site1::AbstractSite, site2::AbstractSite)
    lat = lattice(curr)
    s1 = resolve_site(lat, site1)
    s2 = resolve_site(lat, site2)
    s1 === nothing && return
    s2 === nothing && return
    i = s1.index
    j = s2.index
    if i == j
        @warn "matching source and destination sites"
        return
    end
    curr.currents[i, j] = rhs
    curr.currents[j, i] = -rhs
end

function Base.getindex(curr::Currents, lvm::LatticeValue{Bool})
    check_samesites(curr, lvm)
    indices = findall(lvm.values)
    Currents(curr.lat[lvm], curr.currents[indices, indices])
end

function Base.iterate(curr::Currents{<:SparseMatrixCSC}, state=(1, findnz(curr.currents)))
    ind, (Is, Js, Vs) = state
    ind > length(Is) && return nothing
    i = Is[ind]
    j = Js[ind]
    if j ≤ i
        return iterate(curr, (ind + 1, (Is, Js, Vs)))
    end
    val = Vs[ind]
    return (curr.lat[Is[ind]] => curr.lat[j], val), (ind + 1, (Is, Js, Vs))
end

Base.summary(io::IO, ::Currents{T}) where T = print(io, "Currents{", T, "}")

for f in (:+, :-)
    @eval function ($f)(curr::Currents, curr2::Currents)
        check_samesites(curr, curr2)
        Currents(lattice(curr), ($f)(curr.currents, curr2.currents))
    end
end
for f in (:*, :/)
    @eval ($f)(curr::Currents, num::Number) = Currents(lattice(curr), ($f)(curr.currents, num))
end
*(num::Number, curr::Currents) = curr * num

"""
    Currents(currents[, adjacency_matrix])

Creates a `Currents` instance for `currents`.

## Arguments:
- `currents`: The `AbstractCurrents` object to be turned into `Currents`. That might be time-consuming,
    because  this requires evaluation of the current between all pairs.
- `adjacency_matrix`: If provided, the current will be evaluated only between adjacent sites.
"""
function Currents(curr::AbstractCurrents)
    lat = lattice(curr)
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    for i in eachindex(lat), j in 1:i-1
        ij_curr = curr[i, j]
        abs(ij_curr) < 1e-10 && continue
        push!(Is, i, j)
        push!(Js, j, i)
        push!(Vs, ij_curr, -ij_curr)
    end
    mat = sparse(Is, Js, Vs, length(lat), length(lat))
    return Currents(lat, mat)
end
function Currents(curr::AbstractCurrents, bonds::AbstractBonds)
    lat = lattice(curr)
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    new_bonds = adapt_bonds(bonds, lat)
    for (s1, s2) in new_bonds
        i = s1.index
        j = s2.index
        ij_curr = curr[i, j]
        abs(ij_curr) < 1e-10 && continue
        push!(Is, i, j)
        push!(Js, j, i)
        push!(Vs, ij_curr, -ij_curr)
    end
    mat = sparse(Is, Js, Vs, length(lat), length(lat))
    return Currents(lat, mat)
end

_reorder(p::Pair, ::Nothing) = p
_reorder(p::Pair, by::Function) = by(p[1]) < by(p[2]) ? p : reverse(p)
_mulorder(::Pair, ::Nothing) = 1
_mulorder(p::Pair, by::Function) = by(p[1]) < by(p[2]) ? 1 : -1

"""
    mapgroup_currents(f, group, currents[; sort=false, sortpairsby])

Find the current between all possible pairs of sites, apply `f` to every site pair and
group the result by value of `f`,

## Arguments:
- `f`: This function will be applied to all site pairs. Must accept two `AbstractSite`s.
- `group`: This function will be used to group the current values for pairs with the same mapped value. Must accept a `Vector` of numbers.
- `currents`: The `AbstractCurrents` object to process.

## Keyword arguments:
- `sortresults`: if true, the output arrays will be sorted by results of `f`.
- `sortpairsby`: if provided, the sites in each pair will be sorted by this function.
    Must accept one `AbstractSite`; by default the order of the sites in the pair matches
    their order in the lattice. The sign of the current will match the site order.
"""
function mapgroup_currents(f::Function, group::Function, curr::AbstractCurrents;
        sortresults::Bool=false, sortpairsby::Nullable{Function}=nothing)
    ms = [f(_reorder(pair, sortpairsby)...) for (pair, _) in curr]
    cs = [val * _mulorder(pair, sortpairsby) for (pair, val) in curr]
    new_ms = unique(ms)
    new_cs = [group(cs[ms .== m]) for m in new_ms]
    if sortresults
        perm = sortperm(new_ms)
        permute!(new_ms, perm)
        permute!(new_cs, perm)
    end
    return new_ms, new_cs
end

@recipe function f(curr::AbstractCurrents; showsites=false)
    lat = lattice(curr)
    dims(lat) != 2 && error("2D lattice expected")
    Pts = SVector{2,Float64}[]
    Vs = Float64[]
    nans = SVector(NaN, NaN)
    for ((site1, site2), val) in curr
        abs(val) < 1e-10 && continue
        if val < 0
            site1, site2 = site2, site1
            val = -val
        end
        v1 = site1.coords
        v2 = site2.coords
        d = normalize(v2 - v1)
        o = SVector(d[2], -d[1])
        push!(Pts, v1, v2, v2 - 0.15d - 0.05o, v2 - 0.15d + 0.05o, v2, nans)
        push!(Vs, val, val, val, val, val, NaN)
    end
    @series begin
        seriestype := :path
        seriescolor --> :tempo
        linewidth --> 2.5
        line_z --> Vs
        [p[1] for p in Pts], [p[2] for p in Pts]
    end
    showsites && @series begin
        seriestype := :scatter
        markersize := 1.5
        markercolor := :gray
        markeralpha := 0.8
        markerstrokewidth := 0
        label := ""
        lattice(curr), :sites
    end
end
