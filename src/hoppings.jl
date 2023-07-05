using StaticArrays
import Base: ==

abstract type Boundary end
function shift_site(bc::Boundary, l::Lattice, site::LatticeSite{N}) where N
    lspan = l.lattice_size[i]
    offset = div(site.unit_cell[bc.i] - 1, lspan, RoundDown)
    offset == 0 && return site, 1
    dv = one_hot(bc.i, N)
    factor = 1.
    for _ in 1:abs(offset)
        if offset > 0
            site = displace_site(l, site, -dv * offset)
            factor *= phase_fshift(bc, l, site)
        else
            factor *= phase_fshift(bc, l, site)'
            site = displace_site(l, site, dv * offset)
        end
    end
    factor, site
end

struct TwistedBoundary <: Boundary
    i::Int
    Θ::Float64
end
PeriodicBoundary(i) = TwistedBoundary(i, 0)
phase_fshift(::Lattice, ::LatticeSite, bc::TwistedBoundary) = exp(im * bc.Θ)

struct FunctionBoundaryCondition{F<:Function} <: Boundary
    i::Int
    condition::F
end
phase(bc::Boundary, ::Lattice, ls::LatticeSite) = bc.condition(ls)

struct BoundaryConditions{CondsTuple}
    conds::CondsTuple
    BoundaryConditions(bcs::BT) where BT<:Tuple{<:Boundary} =
        new{BT}(bcs)
end
BoundaryConditions(conds...) = BoundaryConditions(conds)

struct SingleBond{N}
    site1::LatticeSite{N}
    site2::LatticeSite{N}
end

"""
    Bonds

A struct representing bonds in some direction in a lattice.
- `site_indices`: a `NTuple{2, Int}` with indices of sites connected by the bond.
- `translate_uc`: the unit cell offset.
"""
struct Bonds
    site_indices::Tuple{Int,Int}
    translate_uc::Vector{Int}
    function Bonds(site_indices, translate_uc, pbc)
        N = size(hop_operator)[1]
        new{N}(site_indices, translate_uc)
    end
end

_tr_and_pbc(tr_vc::Vector, pbc::Vector) = (tr_vc, pbc)
_tr_and_pbc(tr_vc::Vector, pbc::Bool) = (tr_vc, fill(pbc, length(tr_vc)))
_tr_and_pbc(pbc::Vector) = (zeros(Int, length(pbc)), pbc)
_tr_and_pbc(pbc::Bool) = ([0], [pbc])
function _tr_and_pbc(axis_no::Int, pbc::Vector)
    newdim = max(axis_no, length(pbc))
    tr_vc = zeros(Int, newdim)
    tr_vc[axis_no] = 1
    newpbc = zeros(Bool, newdim)
    newpbc[eachindex(pbc)] = pbc
    return (tr_vc, newpbc)
end
function _tr_and_pbc(axis_no::Int, pbc::Bool)
    tr_vc = zeros(Int, axis_no)
    tr_vc[axis_no] = 1
    return (tr_vc, fill(pbc, axis_no))
end
_get_site_indices(a) = error("Cannot convert object of type $(typeof(a)) to bond indices")
_get_site_indices(i::Int) = (i, i)
_get_site_indices(t::NTuple{2, Int}) = t

"""
    DirectedBonds(kwargs...)

A convenient constructor for a `DirectedBonds` object. `hop_operator` can be either a matrix or a number
(in that case a 1×1 matrix will be created automatically)

**Keyword arguments:**
- `site_indices`: a `NTuple{2, Int}` (or `Int` if they are equal) with indices of sites connected by the bond. `(1, 1)` by default.
- `translate_uc`: the unit cell offset. Zeros by default.
- `axis`: overrides `translate_uc` and sets its components to zero on all axes except given.
- `pbc`: a vector of boolean values indicating if the bonds should be applied periodically over each axis.
Can also be a single boolean, which will set all elements of the vector to given value. `false` by default.

If `site_indices` are equal and `translate_uc` is zero, this means that the bond connects each site with itself,
in which case an error will be thrown.
Note that the dimension count for the hopping is dynamic and will automatically change during runtime.
"""
function Bonds(hop_operator=1; site_indices=1, pbc=false, kw...)
    site_indices = _get_site_indices(site_indices)::NTuple{2, Int}
    if :axis in keys(kw)
        tr_vc, pbc = _tr_and_pbc(kw[:axis], pbc)
    elseif :translate_uc in keys(kw)
        tr_vc, pbc = _tr_and_pbc(kw[:translate_uc], pbc)
    elseif site_indices[1] != site_indices[2]
        tr_vc, pbc = _tr_and_pbc(pbc)
    end
    if iszero(tr_vc) && ==(site_indices...)
        throw(ArgumentError("hopping connects site to itself"))
    end
    Bonds(site_indices, tr_vc, pbc)
end

==(h1::Bonds, h2::Bonds) =
    all(getproperty(h1, fn) == getproperty(h2, fn) for fn in fieldnames(Bonds))

function Base.show(io::IO, m::MIME"text/plain", hop::Bonds)
    println(io, "Bonds connect site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.translate_uc)")
    if !iszero(hop.translate_uc)
        print(io, "Boundary conditions: ")
        if all(hop.pbc .| (hop.translate_uc .== 0))
            println(io, "periodic")
        else
            p_axes = [i for i in eachindex(hop.translate_uc) if hop.pbc[i] && hop.translate_uc[i] != 0]
            if length(p_axes) == 0
                println(io, "open")
            elseif length(p_axes) == 1
                println(io, "periodic at axis $(only(p_axes))")
            else
                println(io, "periodic at axes $(join(p_axes, "×"))")
            end
        end
    end
end
dims(h::Bonds) = length(h.translate_uc)
dims_internal(::Bonds{N}) where N = N

"""
    radius_vector(l::Lattice, hop::Hopping)
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macro cell to minimize the distance between them).
"""
function radius_vector(l::Lattice, hop::Bonds)
    i, j = hop.site_indices
    bravais(l).basis[:, j] - bravais(l).basis[:, i] + bravais(l).translation_vectors * hop.translate_uc
end

"""
    promote_dims!(h::Hopping, ndims::Int)

Changes dimension count of `hopping` to `ndims` if possible.
"""
function promote_dims!(h::Bonds, ndims::Int)
    if ndims ≥ dims(h)
        append!(h.pbc, fill(false, ndims - dims(h)))
        append!(h.translate_uc, fill(0, ndims - dims(h)))
    else
        for _ in 1:dims(h)-ndims
            if h.translate_uc[end] == 0
                pop!(h.translate_uc)
                pop!(h.pbc)
            else
                throw(ArgumentError("cannot shrink hopping dims, non-zero translation found"))
            end
        end
    end
    h
end

"""
    hopping_dest(l::Lattice, hop::Hopping, site::LatticeSite)

Finds the destination site of the `hop` hopping, given the lattice and the source site `site`.

Returns a tuple containing the destination site and a `SVector` with integer numbers describing the macrocell shift.
"""
function hopping_dest(l::Lattice, hop::Bonds, site::LatticeSite)
    site.basis_index != hop.site_indices[1] && return nothing
    new_uc = site.unit_cell + hop.translate_uc
    resid = fld.(new_uc .- 1, size(l))
    all(@. hop.pbc | (resid == 0)) || return nothing
    LatticeSite(mod.(new_uc .- 1, l.lattice_size) .+ 1, hop.site_indices[2], l), resid
end

check_lattice_fits(::Any, ::Lattice) = nothing
@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
function hoppings(selector, lb::LatticeBasis, hops::Bonds...; field::AbstractField=NoField())
    l = lb.latt
    check_lattice_fits(selector, l)
    promote_dims!(hop, dims(l))
    trv = radius_vector(l, hop)
    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    for (i, site1) in enumerate(l)
        for hop in hops
            dst = hopping_dest(l, hop, site1)
            dst === nothing && continue
            site2, _ = dst
            j = site_index(l, site2)
            j === nothing && continue
            !_get_bool_value(selector, l, site1, site2) && continue

            p1 = site1.coords
            pmod = exp(-2π * im * path_integral(field, p1, p1 + trv))
            !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
            push!(Is, i, j)
            push!(Js, j, i)
            push!(Vs, pmod, pmod')
        end
    end
    sparse(Is, Js, Vs, length(lb), length(lb))
end

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
hoppings(lb::LatticeBasis, hops::Bonds...; field::AbstractField=NoField()) =
    hoppings(nothing, lb, hops...; field=field)

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
struct DomainsSelector <: AbstractGraph
    domains::LatticeValue
end
lattice(ps::DomainsSelector) = lattice(ps.domains)
match(ps::DomainsSelector, site1::LatticeSite, site2::LatticeSite) =
    ps.domains[site1] == ps.domains[site2]

"""
    PairLhsSelector(lhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `lhs` is true on the first site of the pair.
"""
struct PairLhsSelector <: AbstractGraph
    lhs::LatticeValue
end
lattice(ps::PairLhsSelector) = lattice(ps.lhs)
match(ps::PairLhsSelector, site1::LatticeSite, ::LatticeSite) = ps.lhs[site1]

"""
    PairRhsSelector(rhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `rhs` is true on the first site of the pair.
"""
struct PairRhsSelector <: AbstractGraph
    rhs::LatticeValue
end
lattice(ps::PairRhsSelector) = lattice(ps.rhs)
match(ps::PairRhsSelector, ::LatticeSite, site2::LatticeSite) = ps.rhs[site2]

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
    bonds(op::LatticeOperator)

Generates a `PairSet` for the provided operator.
"""
function bonds(op::LatticeOperator)
    matrix = Bool[!iszero(op[i, j])
                  for i in 1:length(lattice(op)), j in 1:length(lattice(op))]
    return PairSet(lattice(op), matrix)
end

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
