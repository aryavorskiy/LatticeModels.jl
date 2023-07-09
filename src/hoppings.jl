using StaticArrays

abstract type Boundary end

struct TwistedBoundary <: Boundary
    i::Int
    Θ::Float64
end
PeriodicBoundary(i) = TwistedBoundary(i, 0)
function shift_site(bc::TwistedBoundary, l::Lattice, site::LatticeSite{N}) where N
    ret = 1., site
    bc.i > dims(l) && return ret
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return ret
    dv = one_hot(bc.i, Val(N))
    site = shift_site(l, site, -dv * offset * lspan)
    exp(im * bc.Θ * offset), site
end

struct FunctionBoundary{F<:Function} <: Boundary
    condition::F
    i::Int
end

shift_site(js::SVector{N}, l::Lattice, site::LatticeSite{N}) where N =
    LatticeSite(site.unit_cell + js, site.basis_index, site.coords + bravais(l).translation_vectors * js)

function shift_site(bc::FunctionBoundary, l::Lattice, site::LatticeSite{N}) where N
    ret = 1., site
    bc.i > dims(l) && return ret
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return ret
    dv = one_hot(bc.i, Val(N))
    factor = 1.
    for _ in 1:abs(offset)
        if offset > 0
            site = shift_site(-dv * lspan, l, site)
            factor *= bc.condition(ls)
        else
            factor *= bc.condition(ls)'
            site = shift_site(dv * lspan, l, site)
        end
    end
    factor, site
end

struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    function BoundaryConditions(bcs::CondsTuple) where CondsTuple<:NTuple{N, <:Boundary} where N
        @assert allunique(bc.i for bc in bcs)
        new{CondsTuple}(bcs)
    end
end
_extract_boundary_conditions(b::Boundary) = b
function _extract_boundary_conditions(pb::Pair{Int, Bool})
    !pb.second && error("")
    PeriodicBoundary(pb.first)
end
_extract_boundary_conditions(pb::Pair{Int, <:Real}) = TwistedBoundary(pb...)
BoundaryConditions(args...) = BoundaryConditions(_extract_boundary_conditions.(args))
function shift_site(bcs::BoundaryConditions, l::Lattice, site::LatticeSite)
    factor = 1.
    for bc in bcs.bcs
        new_factor, site = shift_site(bc, l, site)
        factor *= new_factor
    end
    factor, site
end

const SingleBond{N} = Pair{LatticeSite{N}, LatticeSite{N}}

"""
    Bonds{N}

A struct representing bonds in some direction in a lattice.
- `site_indices`: a `NTuple{2, Int}` with indices of sites connected by the bond.
- `translate_uc`: the unit cell offset.
"""
struct Bonds{N}
    site_indices::Pair{Int,Int}
    translate_uc::SVector{N, Int}
end

"""
    Bonds(site_indices; kwargs...)

A convenient constructor for a `Bonds` object.
`site_indices` is a `::Int => ::Int` pair with indices of sites connected by the bond; `1 => 1` is the default value.

**Keyword arguments:**
- `translate_uc`: the unit cell offset, a `Tuple` or a `SVector`. Zeros by default.
- `axis`: overrides `translate_uc` and sets its components to zero on all axes except given.

If `site_indices` are equal and `translate_uc` is zero, this means that the bond connects each site with itself,
in which case an error will be thrown.
Note that the dimension count for the bond is static, it is automatically compatible to higher-dimensional lattices.
"""
function Bonds(site_indices=Pair(1, 1); kw...)
    tr_vc = if :axis in keys(kw)
        one_hot(kw[:axis], kw[:axis])
    elseif :translate_uc in keys(kw)
        kw[:translate_uc]
    elseif site_indices[1] != site_indices[2]
        SA[0]
    end
    if iszero(tr_vc) && ==(site_indices...)
        throw(ArgumentError("bond connects site to itself"))
    end
    Bonds(site_indices, tr_vc)
end

Base.:(==)(h1::Bonds, h2::Bonds) =
    all(getproperty(h1, fn) == getproperty(h2, fn) for fn in fieldnames(Bonds))

function Base.show(io::IO, ::MIME"text/plain", hop::Bonds)
    println(io, "Bonds connect site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.translate_uc)")
end
dims(::Bonds{N}) where N = N

"""
    radius_vector(l::Lattice, hop::Hopping)
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macro cell to minimize the distance between them).
"""
function radius_vector(l::Lattice, hop::Bonds)
    i, j = hop.site_indices
    bravais(l).basis[:, j] - bravais(l).basis[:, i] +
     mm_assuming_zeros(bravais(l).translation_vectors, hop.translate_uc)
end

check_lattice_fits(::Any, ::Lattice) = nothing
@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
# @inline _get_bool_value(g::AbstractGraph, ::Lattice, site1::LatticeSite, site2::LatticeSite) =
#     match(g, site1, site2)

function add_hoppings!(builder, selector, l::Lattice, op, bond::Bonds,
        field::AbstractField, boundaries::BoundaryConditions)
    dims(bond) > dims(l) && error("Incompatible dims")
    trv = radius_vector(l, bond)
    for site1 in l
        site1.basis_index != hop.site_indices[1] && continue
        site2 = LatticeSite(add_assuming_zeros(site1.unit_cell, hop.translate_uc),
            hop.site_indices[2], site1.coords + trv)
        add_hoppings!(builder, selector, l, op, site1 => site2, field, boundaries)
    end
end

function add_hoppings!(builder, selector, l::Lattice, bond::SingleBond,
        field::AbstractField, boundaries::BoundaryConditions)
    site1, site2 = bond
    p1 = site1.coords
    p2 = site2.coords
    factor, site2 = shift_site(boundaries, l, site2)
    i = @inline site_index(l, site1)
    j = @inline site_index(l, site2)
    i === nothing && return
    j === nothing && return
    !_get_bool_value(selector, l, site1, site2) && return
    total_factor = exp(-2π * im * line_integral(field, p1, p2)) * factor
    !isfinite(total_factor) && error("got NaN or Inf when finding the phase factor")
    increment!(builder, total_factor, i, j)
    increment!(builder, total_factor', j, i)
end

function hoppings(selector, lb::LatticeBasis, hops...;
        field::AbstractField=NoField(), boundaries::BoundaryConditions=BoundaryConditions())
    l = lb.latt
    check_lattice_fits(selector, l)
    builder = SparseMatrixBuilder((length(lb), length(lb)))
    for hop in hops
        add_hoppings!(builder, selector, l, hop, field, boundaries)
    end
    Operator(lb, to_matrix(builder))
end
hoppings(selector, lb::LatticeBasis; kw...) =
    hoppings(selector, lb, default_bonds(lb.latt)...; kw...)
hoppings(selector, l::Lattice, args...; kw...) =
    hoppings(selector, LatticeBasis(l), args...; kw...)

@doc raw"""
    hoppings([f, ]lattice::Lattice, hopping::Hopping[, field::AbstractField])

Creates a hopping operator:
$$\hat{A} = \sum_{pairs} \hat{c}^\dagger_j \hat{c}_i + h. c.$$

Arguments:
- `f`: a function that takes a `Lattice` and two `LatticeSite`s, returns whether this pair should be included.
Can also be a `PairSelector`,
- `lattice`: the lattice to create the operator on.
- `hopping`: the `Hopping` object describing the site pairs and the $\hat{t}$ operator.
- `field`: the `AbstractField` object that defines the magnetic field to generate phase factors using Peierls substitution.
"""
hoppings(l, hops::Bonds...; field::AbstractField=NoField(), boundaries=BoundaryConditions()) =
    hoppings(nothing, l, hops...; field=field, boundaries=boundaries)

"""
    AbstractGraph <: Function

A function-like object that accepts a `Lattice` and two `LatticeSite`s and returns whether the site pair is *selected* or not.
Implements a built-in sanity check algorithm to make sure the pair set was defined on a correct lattice.

Define these functions for all subclasses:
- `LatticeModels.lattice(::YourSelector)` must return the lattice your selector was defined on.
- `LatticeModels.match(::YourSelector, site1::LatticeSite, site2::LatticeSite)` must return whether the `(site1, site2)` pair is *selected*.
"""
abstract type AbstractGraph <: Function end
match(ps::AbstractGraph, ::LatticeSite, ::LatticeSite) =
    error("match(::$(typeof(ps)), ::LatticeSite, ::LatticeSite) must be explicitly implemented")
lattice(ps::AbstractGraph) = error("lattice(::$(typeof(ps))) must be explicitly implemented")
check_lattice_fits(ps::AbstractGraph, l::Lattice) = check_is_sublattice(l, lattice(ps))
(ps::AbstractGraph)(::Lattice, site1::LatticeSite, site2::LatticeSite) = match(ps, site1, site2)

"""
    Domains(domains::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `domains` is the same on two sites.
"""
struct Domains{LT} <: AbstractGraph
    domains::LT
    Domains(domains::LT) where LT<:LatticeValue = new{LT}(domains)
end
lattice(ps::Domains) = lattice(ps.domains)
match(ps::Domains, site1::LatticeSite, site2::LatticeSite) =
    ps.domains[site1] == ps.domains[site2]

"""
    PairLhsGraph(lhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `lhs` is true on the first site of the pair.
"""
struct PairLhsGraph <: AbstractGraph
    lhs::LatticeValue
end
lattice(ps::PairLhsGraph) = lattice(ps.lhs)
match(ps::PairLhsGraph, site1::LatticeSite, ::LatticeSite) = ps.lhs[site1]

"""
    PairRhsGraph(rhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `rhs` is true on the first site of the pair.
"""
struct PairRhsGraph <: AbstractGraph
    rhs::LatticeValue
end
lattice(ps::PairRhsGraph) = lattice(ps.rhs)
match(ps::PairRhsGraph, ::LatticeSite, site2::LatticeSite) = ps.rhs[site2]

struct InvertedGraph{GT} <: AbstractGraph
    graph::GT
end
lattice(ig::InvertedGraph) = lattice(ig.graph)
match(ig::InvertedGraph, site1, site2) = match(ig.graph, site1, site2)
Base.:(!)(g::AbstractGraph) = InvertedGraph(g)
Base.:(!)(g::InvertedGraph) = g.graph
