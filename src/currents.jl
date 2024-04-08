import Base: getindex, +, -, *, /
using LinearAlgebra

const CURRENTS_EPS = 1e-10

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
    return (lat[i] => lat[j], curr[i, j]), i => j
end

"""
    SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents

A lazy wrapper for a `Currents` object representing the same currents but on a smaller lattice.
"""
struct SubCurrents{CT} <: AbstractCurrents
    parent_currents::CT
    lat::AbstractLattice
    indices::Vector{Int}
    function SubCurrents(parent_currents::CT, indices::Vector{Int}) where {CT<:AbstractCurrents}
        lat = lattice(parent_currents)
        new{CT}(parent_currents, lat[indices], indices)
    end
end

Base.getindex(scurr::SubCurrents, i::Int, j::Int) = scurr.parent_currents[scurr.indices[i], scurr.indices[j]]
lattice(scurr::SubCurrents) = scurr.lat
function Base.getindex(curr::AbstractCurrents, any)
    inds = to_inds(lattice(curr), any)
    SubCurrents(curr, inds)
end

function Base.summary(io::IO, scurr::SubCurrents)
    print(io, "SubCurrents of ")
    summary(io, scurr.parent_currents)
end

"""
    currentsfrom(currents, src)

Create a `LatticeValue` object with the currents from `src` region to all other sites.

## Arguments
- `currents`: The `AbstractCurrents` object to process.
- `src`: The source region. Can be a site/collection of sites or a `LatticeValue{Bool}` mask.
"""
function currentsfrom(curr::AbstractCurrents, src)
    lat = lattice(curr)
    is = to_inds(lat, src)
    LatticeValue(lat, Float64[j in is ? 0 : sum(curr[i, j] for i in is) for j in eachindex(lat)])
end

"""
    currentsfromto(currents, src[, dst])

Finds the total current from `src` to `dst` regions. If `dst` is not provided, the current
from `src` to all other sites is returned.

## Arguments
- `currents`: The `AbstractCurrents` object to process.
- `src`: The source region.
- `dst`: The destination region.

Both `src` and `dst` can be a site/collection of sites or a `LatticeValue{Bool}` mask.
"""
function currentsfromto(curr::AbstractCurrents, src, dst=nothing)
    lat = lattice(curr)
    src_inds = to_inds(lat, src)
    dst_inds = dst === nothing ? setdiff(eachindex(lat), src_inds) : to_inds(lat, dst)
    sum(curr[i, j] for i in src_inds, j in dst_inds)
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

function Base.getindex(curr::Currents, any)
    inds = to_inds(lattice(curr), any)
    Currents(curr.lat[inds], curr.currents[inds, inds])
end

function SparseArrays.findnz(curr::AbstractCurrents)
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    lat = lattice(curr)
    for j in 1:length(lat), i in 1:j-1
        val = curr[i, j]
        abs(val) < CURRENTS_EPS && continue
        push!(Is, i)
        push!(Js, j)
        push!(Vs, val)
    end
    return Is, Js, Vs
end
function SparseArrays.findnz(curr::Currents{<:SparseMatrixCSC})
    Is, Js, Vs = findnz(curr.currents)
    mask = Is .< Js
    return Is[mask], Js[mask], Vs[mask]
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

## Arguments
- `currents`: The `AbstractCurrents` object to be turned into `Currents`. That might be time-consuming,
    because  this requires evaluation of the current between all pairs.
- `adjacency_matrix`: If provided, the current will be evaluated only between adjacent sites.

## Examples
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(4, 4); site1, site2 = lat[1:2];

julia> H0 = tightbinding_hamiltonian(lat); psi = groundstate(H0);

julia> H1 = tightbinding_hamiltonian(lat, field=LandauGauge(0.1));

julia> currents = DensityCurrents(H1, psi)
Density currents for system:
One particle on 16-site SquareLattice in 2D space

julia> c2 = Currents(currents)
Currents{SparseArrays.SparseMatrixCSC{Float64, Int64}} on 16-site SquareLattice in 2D space

julia> c2[site1, site2] ≈ currents[site1, site2]
true
```
"""
function Currents(curr::AbstractCurrents)
    lat = lattice(curr)
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    for i in eachindex(lat), j in 1:i-1
        ij_curr = curr[i, j]
        abs(ij_curr) < CURRENTS_EPS && continue
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
        abs(ij_curr) < CURRENTS_EPS && continue
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
    mapgroup_currents(f, group, currents[; sortresults=false, sortpairsby])

Find the current between all possible pairs of sites, apply `f` to every site pair and
group the result by value of `f`,

## Arguments
- `f`: This function will be applied to all site pairs. Must accept two `AbstractSite`s.
- `group`: This function will be used to group the current values for pairs with the same mapped value. Must accept a `Vector` of numbers.
- `currents`: The `AbstractCurrents` object to process.

## Keyword arguments
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
