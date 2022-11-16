import Base: getindex
using LinearAlgebra

"""
    AbstractCurrents

Supertype for all type representing currents-like values on a lattice.
Subtypes must implement `current_lambda` and `lattice` functions.
"""
abstract type AbstractCurrents end

"""
    current_lambda(::AbstractCurrents)

Returns a function that takes two integer indices of sites in a lattice and returns the current between these two sites.
"""
current_lambda(::T) where {T<:AbstractCurrents} = error("current_lambda(::$T) must be explicitly implemented")

"""
    lattice(::AbstractCurrents)

Gets the lattice where the given `AbstractCurrents` object is defined.
"""
lattice(::T) where {T<:AbstractCurrents} = error("lattice(::$T) must be explicitly implemented")

"""
    DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents <: AbstractCurrents
    hamiltonian::LatticeOperator
    density::LatticeOperator

    """
        DensityCurrents(hamiltonian, density_mat)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `density_mat`.
    """
    function DensityCurrents(ham::LatticeOperator, dens::LatticeOperator)
        ham.basis != dens.basis && error("basis mismatch")
        new(ham, dens)
    end
end

current_lambda(curr::DensityCurrents) =
    (i::Int, j::Int) -> 2imag(tr(curr.density[i, j] * curr.hamiltonian[j, i]))
lattice(curr::DensityCurrents) = curr.hamiltonian.basis.lattice

"""
    SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents

A lazy wrapper for a `SubCurrents` object that representing the same currents but on a smaller lattice.
"""
struct SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents
    parent_currents::CT
    lattice::Lattice
    indices::Vector{Int}
    function SubCurrents(parent_currents::CT, indices::Vector{Int}) where {CT}
        l = lattice(parent_currents)
        m = zeros(Bool, length(l))
        m[indices] .= true
        new_mask = zero(l.mask)
        new_mask[l.mask] = m
        new{CT}(parent_currents,
            Lattice(lattice_type(l), size(l), bravais(l), new_mask), indices)
    end
end

function current_lambda(scurr::SubCurrents)
    in_fn = current_lambda(scurr.parent_currents)
    (i::Int, j::Int) -> in_fn(scurr.indices[i], scurr.indices[j])
end
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

lattice(mcurr::MaterializedCurrents) = mcurr.lattice
current_lambda(mcurr::MaterializedCurrents) = (i::Int, j::Int) -> mcurr.currents[i, j]

function getindex(curr::AbstractCurrents, lvm::LatticeValue{Bool})
    lattice(curr) != lattice(lvm) && error("lattice mismatch")
    indices = findall(lvm.values)
    SubCurrents(curr, indices)
end

function getindex(curr::MaterializedCurrents, lvm::LatticeValue{Bool})
    lattice(curr) != lattice(lvm) && error("lattice mismatch")
    indices = findall(lvm.values)
    MaterializedCurrents(curr.lattice[lvm], curr.currents[indices, indices])
end

"""
    materialize([function, ]currents)

Creates a `MaterializedCurrents` instance for `currents`.

If `function` is provided, it must accept a `Lattice` and two `LatticeSite`s and return if the current between this site must be calculated or not.
This can be useful to avoid exsessive calculations.
"""
function materialize(curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    curr_fn = current_lambda(curr)
    for i in 1:length(l), j in 1:i-1
        ij_curr = curr_fn(i, j)
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end

function materialize(f::Function, curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    curr_fn = current_lambda(curr)
    for i in 1:length(l), j in 1:i-1
        if !f(l, i, j)
            continue
        end
        ij_curr = curr_fn(i, j)
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end

"""
    pairs_by_adjacent(bonds)

A selector function used for hopping operator definition or currents materialization.

Takes a `BondSet` and generates a lambda which accepts a lattice and two integer site indices,
returning whether the two sites are connected by the `bonds`.
"""
pairs_by_adjacent(bs::BondSet) =
    (l::Lattice, i::Int, j::Int) -> is_adjacent(bs, l[i], l[j])

"""
    pairs_by_distance(function)

A selector function used for hopping operator definition or currents materialization.

Takes a function and generates a lambda which accepts a lattice and two integer site indices,
returning whether `function` applied to distence between the two sites returned `true`.
"""
pairs_by_distance(f) =
    (l::Lattice, i::Int, j::Int) -> f(norm(radius_vector(l, l[i], l[j])))::Bool

"""
    map_currents(map_fn, currents[; aggr_fn, sort])

Accepts a function that takes a `Lattice` and two `LatticeSite`s and returns any value.
Applies `map_fn` to every site pair and returns two `Vector`s: one with currents, one with results of `map_fn`.

**Keyword arguments:**
- `aggr_fn`: if a function is provided, all currents with the same mapped value will be aggregated into one value.
For example, if `aggr_fn=(x -> mean(abs.(x)))`, and `map_fn` finds the distance between the sites,
the returned lists will store the distance between sites and the average absolute current between sites with such distance.
"""
function map_currents(f::Function, curr::AbstractCurrents; aggr_fn::Union{Nothing, Function}=nothing, sorted::Bool=false)
    l = lattice(curr)
    curr_fn = current_lambda(curr)
    cs = Float64[]
    ms = only(Base.return_types(f, (Lattice, LatticeSite, LatticeSite)))[]
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            if site1 > site2
                push!(cs, curr_fn(i, j))
                push!(ms, f(l, site1, site2))
            else
                break
            end
            j += 1
        end
        i += 1
    end
    if aggr_fn !== nothing
        ms_set = collect(Set(ms))
        cs = [aggr_fn(cs[ms .== m]) for m in ms_set]
        ms = ms_set
    end
    if sorted
        perm = sortperm(ms)
        return ms[perm], cs[perm]
    end
    ms, cs
end

@recipe function f(curr::AbstractCurrents)
    l = lattice(curr)
    dims(l) != 2 && error("2D lattice expected")
    Xs = Float64[]
    Ys = Float64[]
    Qs = NTuple{2,Float64}[]
    curr_fn = current_lambda(curr)
    arrows_scale --> 1
    arrows_rtol --> 1e-2
    seriestype := :quiver
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            j ≥ i && continue
            ij_curr = curr_fn(i, j)::Real
            crd = site_coords(l, (ij_curr > 0 ? site1 : site2))
            vc = radius_vector(l, site2, site1)
            vc_n = norm(vc)
            if vc_n < abs(ij_curr * plotattributes[:arrows_scale] / plotattributes[:arrows_rtol])
                push!(Xs, crd[1])
                push!(Ys, crd[2])
                push!(Qs, Tuple(vc * (ij_curr * plotattributes[:arrows_scale] / vc_n)))
            end
            j += 1
        end
        i += 1
    end
    quiver := Qs
    Xs, Ys
end
