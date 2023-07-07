using StaticArrays
import Base: ==

abstract type Boundary end

struct TwistedBoundary <: Boundary
    i::Int
    Θ::Float64
end
PeriodicBoundary(i) = TwistedBoundary(i, 0)
function shift_site(bc::TwistedBoundary, l::Lattice, site::LatticeSite{N}) where N
    bc.i > dims(l) && return 1., site
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return 1., site
    dv = one_hot(bc.i, Val(N))
    site = displace_site(l, site, -dv * offset * lspan)
    exp(im * bc.Θ * offset), site
end

struct FunctionBoundary{F<:Function} <: Boundary
    i::Int
    condition::F
end
function shift_site(bc::FunctionBoundary, l::Lattice, site::LatticeSite{N}) where N
    bc.i > dims(l) && return 1., site
    lspan = l.lattice_size[bc.i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return site, 1
    dv = one_hot(bc.i, Val(N))
    factor = 1.
    for _ in 1:abs(offset)
        if offset > 0
            site = displace_site(l, site, -dv * lspan)
            factor *= bc.condition(ls)
        else
            factor *= bc.condition(ls)'
            site = displace_site(l, site, dv * lspan)
        end
    end
    factor, site
end

struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    function BoundaryConditions(bcs::BT) where BT<:NTuple{N, <:Boundary} where N
        @assert allunique(bc.i for bc in bcs)
        new{BT}(bcs)
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

struct SingleBond{N}
    site1::LatticeSite{N}
    site2::LatticeSite{N}
end

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

==(h1::Bonds, h2::Bonds) =
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

"""
    hopping_dest(l::Lattice, hop::Hopping, site::LatticeSite)

Finds the destination site of the `hop` hopping, given the lattice and the source site `site`.

Returns a tuple containing the destination site and a `SVector` with integer numbers describing the macrocell shift.
"""
function hopping_dest(l::Lattice, hop::Bonds, site::LatticeSite)
    site.basis_index != hop.site_indices[1] && return nothing
    new_uc = add_assuming_zeros(site.unit_cell, hop.translate_uc)
    get_site(l, new_uc, hop.site_indices[2])
end

check_lattice_fits(::Any, ::Lattice) = nothing
@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
function hoppings(selector, lb::LatticeBasis, hops::Bonds...; field::AbstractField=NoField(), boundaries::BoundaryConditions=BoundaryConditions())
    l = lb.latt
    check_lattice_fits(selector, l)
    trvs = [radius_vector(l, hop) for hop in hops]
    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    for hop in hops
        dims(hop) > dims(l) && error("Incompatible dims")
    end
    for (i, site1) in enumerate(l)
        for (hi, hop) in enumerate(hops)
            site1.basis_index != hop.site_indices[1] && continue
            site2 = LatticeSite(add_assuming_zeros(site1.unit_cell, hop.translate_uc),
                hop.site_indices[2], site1.coords + trvs[hi])
            factor, site2 = shift_site(boundaries, l, site2)
            j = site_index(l, site2)
            j === nothing && continue
            !_get_bool_value(selector, l, site1, site2) && continue

            p1 = site1.coords
            pmod = exp(-2π * im * line_integral(field, p1, p1 + trvs[hi])) * factor
            !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
            push!(Is, i, j)
            push!(Js, j, i)
            push!(Vs, pmod, pmod')
        end
    end
    mat = sparse(Is, Js, Vs, length(lb), length(lb))
    Operator(lb, mat)
end
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
    DomainsSelector(domains::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `domains` is the same on two sites.
"""
struct DomainsSelector{LT} <: AbstractGraph
    domains::LT
    DomainsSelector(domains::LT) where LT<:LatticeValue = new{LT}(domains)
end
lattice(ps::DomainsSelector) = lattice(ps.domains)
match(ps::DomainsSelector, site1::LatticeSite, site2::LatticeSite) =
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

macro hopping_operator(for_loop::Expr)
    if for_loop.head !== :for
        throw(ArgumentError("expression must be a for loop, not $(for_loop.head)"))
    elseif length(for_loop.args) != 2
        throw(ArgumentError("malformed for loop")) # This should never happen, but still...
    end
    itr::Expr, body::Expr = for_loop.args
    if !Meta.isexpr(itr, :(=), 2)
        throw(ArgumentError("invalid loop iteration specification"))
    end
    itr_vars, lattice_var = itr.args
    if !Meta.isexpr(itr_vars, :tuple, 2)
        throw(ArgumentError("invalid loop iteration variable; must be a LatticeSite 2-tuple"))
    end
    site1_var, site2_var = itr_vars.args
    while body.args[end] isa LineNumberNode
        pop!(body.args)
    end
    dump(body)
    quote
        l = $(esc(lattice_var))
        local matrix = nothing
        for (i, $(esc(site1_var))) in enumerate(l)
            for (j, $(esc(site2_var))) in enumerate(l)
                if i ≥ j
                    continue
                end
                block_res = $(esc(body))
                if block_res !== nothing
                    if matrix === nothing
                        matrix = zero_on_basis(basis(l, block_res))
                    end
                    matrix[i, j] .= block_res
                    matrix[j, i] .= block_res'
                end
            end
        end
        matrix
    end
end

"""
    PairSet{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`PairSet`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `PairSet` which connects sites that were connected by `≤n` bonds of the previous `PairSet`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct PairSet{LT<:Lattice} <: AbstractGraph
    lattice::LT
    bmat::Matrix{Bool}
    function PairSet(l::LT, bmat) where {LT<:Lattice}
        !all(size(bmat) .== length(l)) && error("inconsistent connectivity matrix size")
        new{LT}(l, bmat .| bmat' .| Matrix(I, length(l), length(l)))
    end
    function PairSet(l::Lattice)
        PairSet(l, Matrix(I, length(l), length(l)))
    end
end

lattice(bs::PairSet) = bs.lattice
match(bs::PairSet, site1::LatticeSite, site2::LatticeSite) =
    bs.bmat[site_index(lattice(bs), site1), site_index(lattice(bs), site2)]

"""
    bonds(l::Lattice, hoppings::Hopping...)

Generates a `PairSet` for a given set of `Hopping`s on a given `Lattice`.
"""
function bonds(l::Lattice, hops::Bonds...)
    bs = PairSet(l)
    for h in hops
        promote_dims!(h, dims(l))
    end
    for i in 1:length(l)
        site1 = bs.lattice[i]
        for hop in hops
            dst = hopping_dest(l, hop, site1)
            dst === nothing && continue
            site2 = dst[1]
            j = site_index(l, site2)
            j === nothing && continue
            bs.bmat[i, j] = bs.bmat[j, i] = true
        end
    end
    bs
end

import Base: !, ^, |

function |(bss::PairSet...)
    !allequal(getproperty.(bss, :lattice)) && error("lattice mismatch")
    PairSet(bss[1].lattice, .|(getproperty.(bss, :bmat)...))
end

^(bs1::PairSet, n::Int) = PairSet(bs1.lattice, bs1.bmat^n .!= 0)

!(bs::PairSet) =
    PairSet(bs.lattice, Matrix(I, length(bs.lattice), length(bs.lattice)) .| .!(bs.bmat))

function Base.show(io::IO, m::MIME"text/plain", bs::PairSet)
    println(io, "PairSet with $(count(bs.bmat)) bonds")
    print(io, "on ")
    show(io, m, bs.lattice)
end
