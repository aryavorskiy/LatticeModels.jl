import Base: getindex, +, -, *, /
using LinearAlgebra

"""
    AbstractCurrents

Supertype for all type representing currents-like values on a lattice.
Subtypes must implement `Base.getindex(::Int, ::Int)` and `lattice` functions.
"""
abstract type AbstractCurrents end

"""
    lattice(::AbstractCurrents)

Gets the lattice where the given `AbstractCurrents` object is defined.
"""
lattice(curr::AbstractCurrents) = error("lattice(::$(typeof(curr))) must be explicitly implemented")

Base.getindex(curr::AbstractCurrents, s1::LatticeSite, s2::LatticeSite) =
    curr[site_index(lattice(curr), s1), site_index(lattice(curr), s2)]

"""
    SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents

A lazy wrapper for a `Currents` object representing the same currents but on a smaller lattice.
"""
struct SubCurrents{CT} <: AbstractCurrents
    parent_currents::CT
    lattice::Lattice
    indices::Vector{Int}
    function SubCurrents(parent_currents::CT, indices::Vector{Int}) where {CT<:AbstractCurrents}
        l = lattice(parent_currents)
        m = zeros(Bool, length(l))
        m[indices] .= true
        new_mask = zero(l.mask)
        new_mask[l.mask] = m
        new{CT}(parent_currents,
            Lattice(lattice_type(l), size(l), bravais(l), new_mask), indices)
    end
end

Base.getindex(scurr::SubCurrents, i::Int, j::Int) = scurr.parent_currents[scurr.indices[i], scurr.indices[j]]
lattice(scurr::SubCurrents) = scurr.lattice

"""
    MaterializedCurrents <: AbstractCurrents

A `AbstractCurrents` instance that stores values for all currents explicitly.
"""
struct MaterializedCurrents <: AbstractCurrents
    lattice::Lattice
    currents::Matrix{Float64}
    function MaterializedCurrents(l::Lattice, curs::Matrix{Float64})
        !all(length(l) .== size(curs)) && error("dimension mismatch")
        new(l, curs)
    end
end

MaterializedCurrents(l::Lattice) =
    MaterializedCurrents(l, zeros(length(l), length(l)))

Base.convert(::Type{MaterializedCurrents}, curr::AbstractCurrents) = materialize(curr)
Base.copy(mc::MaterializedCurrents) = MaterializedCurrents(lattice(mc), copy(mc.currents))
Base.zero(mc::MaterializedCurrents) = MaterializedCurrents(lattice(mc), zero(mc.currents))
lattice(mcurr::MaterializedCurrents) = mcurr.lattice

Base.getindex(mcurr::MaterializedCurrents, i::Int, j::Int) = mcurr.currents[i, j]
function Base.getindex(curr::AbstractCurrents, lvm::LatticeValue{Bool})
    check_lattice_match(curr, lvm)
    indices = findall(lvm.values)
    SubCurrents(curr, indices)
end

function Base.getindex(curr::MaterializedCurrents, lvm::LatticeValue{Bool})
    check_lattice_match(curr, lvm)
    indices = findall(lvm.values)
    MaterializedCurrents(curr.lattice[lvm], curr.currents[indices, indices])
end

function _site_indices(l::Lattice, l2::Lattice)
    check_is_sublattice(l, l2)
    [site_index(l, site) for site in l2]
end
_site_indices(l::Lattice, site::LatticeSite) = (site_index(l, site),)
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

for f in (:+, :-)
    @eval function ($f)(curr::MaterializedCurrents, curr2::MaterializedCurrents)
        check_lattice_match(curr, curr2)
        MaterializedCurrents(lattice(curr), ($f)(curr.currents, curr2.currents))
    end
end
for f in (:*, :/)
    @eval ($f)(curr::MaterializedCurrents, num::Number) = MaterializedCurrents(lattice(curr), ($f)(curr.currents, num))
end
*(num::Number, curr::MaterializedCurrents) = curr * num

"""
    materialize([function, ]currents)

Creates a `MaterializedCurrents` instance for `currents`.

If `function` is provided, it must accept a `Lattice` and two `LatticeSite`s and return if the current between this site must be calculated or not.
This can be useful to avoid exsessive calculations.
"""
function materialize(curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    for i in 1:length(l), j in 1:i-1
        ij_curr = curr[i, j]
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end

function materialize(f::Function, curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    for i in 1:length(l), j in 1:i-1
        !f(l, l[i], l[j]) && continue
        ij_curr = curr[i, j]
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end

"""
    pairs_by_distance(f)

A selector function used for hopping operator definition or currents materialization.

Takes a function and generates a lambda which accepts a lattice and two `LatticeSite`s,
returning whether `f` applied to distance between the two sites returned `true`.
"""
pairs_by_distance(f) =
    (l::Lattice, site1::LatticeSite, site2::LatticeSite) ->
        f(norm(radius_vector(l, site1, site2)))::Bool

"""
    map_currents(map_fn, currs::AnstractCurrents[; reduce_fn, sort=false])

Accepts a function that takes a `Lattice` and two `LatticeSite`s and returns any value.
Applies `map_fn` to every site pair and returns two `Vector`s: one with currents, one with results of `map_fn`.

**Keyword arguments:**
- `reduce_fn`: if a function is provided, all currents with the same mapped value will be reduced into one value.
For example, if `aggr_fn=(x -> mean(abs.(x)))`, and `map_fn` finds the distance between the sites,
the returned lists will store the distance between sites and the average absolute current between sites with such distance.
- `sort`: if true, the output arrays will be sorted by mapped value.
"""
function map_currents(f::Function, curr::AbstractCurrents; reduce_fn::Nullable{Function}=nothing, sort::Bool=false)
    l = lattice(curr)
    cs = Float64[]
    ms = only(Base.return_types(f, (Lattice, LatticeSite, LatticeSite)))[]
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            if site1 > site2
                push!(cs, curr[i, j])
                push!(ms, f(l, site1, site2))
            else
                break
            end
        end
    end
    if reduce_fn !== nothing
        ms_set = unique(ms)
        cs = [reduce_fn(cs[ms .== m]) for m in ms_set]
        ms = ms_set
    end
    if sort
        perm = sortperm(ms)
        permute!(ms, perm)
        permute!(cs, perm)
    end
    if eltype(cs) <: Vector
        ms, transpose(hcat(cs...))
    else
        ms, cs
    end
end
