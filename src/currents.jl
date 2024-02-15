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
        l = lattice(curr)
        return (l[1] => l[2], curr[1, 2]), 1 => 2
    end
end
@inline function Base.iterate(curr::AbstractCurrents, pair::Pair{Int, Int})
    l = lattice(curr)
    i, j = pair
    j += 1
    if j > length(l)
        i += 1
        j = i + 1
        j > length(l) && return nothing
    end
    iszerocurrent(curr, i, j) && return iterate(curr, i => j)
    return (l[i] => l[j], curr[i, j]), i => j
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
        l = lattice(parent_currents)
        new{CT}(parent_currents, l[indices], indices)
    end
end

Base.getindex(scurr::SubCurrents, i::Int, j::Int) = scurr.parent_currents[scurr.indices[i], scurr.indices[j]]
lattice(scurr::SubCurrents) = scurr.lattice
function Base.getindex(curr::AbstractCurrents, lvm::LatticeValue{Bool})
    check_samesites(curr, lvm)
    indices = findall(lvm.values)
    SubCurrents(curr, indices)
end

function _site_indices(l::AbstractLattice, l2::AbstractLattice)
    check_issublattice(l2, l)
    return [site_index(l, site) for site in l2]
end
_site_indices(l::AbstractLattice, site::AbstractSite) = (site_index(l, site),)
function currents_from(curr::AbstractCurrents, src)
    l = lattice(curr)
    is = _site_indices(l, src)
    LatticeValue(l, Float64[j in is ? 0 : sum(curr[i, j] for i in is) for j in eachindex(l)])
end
function currents_from_to(curr::AbstractCurrents, src, dst=nothing)
    l = lattice(curr)
    is = _site_indices(l, src)
    js = dst === nothing ? setdiff(eachindex(l), is) : _site_indices(l, dst)
    sum(curr[i, j] for i in is, j in js)
end

"""
    Currents <: AbstractCurrents

A `AbstractCurrents` instance that stores values for all currents explicitly.
"""
struct Currents{T, LT} <: AbstractCurrents
    lattice::LT
    currents::Matrix{Float64}
    function Currents(l::LT, curs::Matrix{T}) where {T, LT}
        @check_size curs :square
        @check_size l size(curs, 1)
        new{T, LT}(l, curs)
    end
end
Currents(l::AbstractLattice) =
    Currents(l, zeros(length(l), length(l)))

Base.convert(::Type{Currents}, curr::AbstractCurrents) = Currents(curr)
Base.copy(mc::Currents) = Currents(lattice(mc), copy(mc.currents))
Base.zero(mc::Currents) = Currents(lattice(mc), zero(mc.currents))
lattice(mcurr::Currents) = mcurr.lattice


Base.:(==)(c1::Currents, c2::Currents) =
    c1.lattice == c2.lattice && c1.currents == c2.currents
Base.isapprox(c1::Currents, c2::Currents; kw...) =
    c1.lattice == c2.lattice && isapprox(c1.currents, c2.currents; kw...)

Base.getindex(mcurr::Currents, i::Int, j::Int) = mcurr.currents[i, j]
function Base.setindex!(curr::Currents, rhs, site1::AbstractSite, site2::AbstractSite)
    l = lattice(curr)
    s1 = resolve_site(l, site1)
    s2 = resolve_site(l, site2)
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
    Currents(curr.lattice[lvm], curr.currents[indices, indices])
end

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
    l = lattice(curr)
    m = Currents(l)
    for i in eachindex(l), j in 1:i-1
        ij_curr = curr[i, j]
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end
function Currents(curr::AbstractCurrents, bonds::AbstractBonds)
    l = lattice(curr)
    m = Currents(l)
    new_bonds = apply_lattice(bonds, l)
    for (s1, s2) in new_bonds
        i = s1.index
        j = s2.index
        ij_curr = curr[i, j]
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
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

@recipe function f(curr::AbstractCurrents)
    l = lattice(curr)
    dims(l) != 2 && error("2-dim lattice expected")
    Xs = Float64[]
    Ys = Float64[]
    Qs = NTuple{2,Float64}[]
    arrows_scale --> 1
    arrows_rtol --> 1e-2
    seriestype := :quiver
    for ((site1, site2), val) in curr
        crd = val > 0 ? site1.coords : site2.coords
        vc = site2.coords - site1.coords
        vc_n = norm(vc)
        if vc_n < abs(val * plotattributes[:arrows_scale] / plotattributes[:arrows_rtol])
            push!(Xs, crd[1])
            push!(Ys, crd[2])
            push!(Qs, Tuple(vc * (val * plotattributes[:arrows_scale] / vc_n)))
        end
    end
    quiver := Qs
    Xs, Ys
end
