using StaticArrays
import Base: copy, show, ==

"""
    Hopping{MT} where {MT<:AbstractMatrix}

A struct representing a bond in a lattice.
- `site_indices`: a `NTuple{2, Int}` with indices of sites connected by the bond.
- `translate_uc`: the unit cell offset.
- `pbc`: a vector of boolean values indicating if the bonds should be applied periodically over each axis.
- `hop_operator`: a matrix of type `MT` representing the operator affecting the internal state.
"""
struct Hopping{MT<:AbstractMatrix}
    site_indices::Tuple{Int,Int}
    translate_uc::Vector{Int}
    pbc::Vector{Bool}
    hop_operator::MT
    function Hopping(site_indices, translate_uc, pbc, hop_operator)
        length(translate_uc) != length(pbc) && error("inconsistent dimensionality")
        new{typeof(hop_operator)}(site_indices, translate_uc, pbc, hop_operator)
    end
end

_wrap_operator(o::Number) = [o;;]
_wrap_operator(o::AbstractMatrix) = o
_tr_and_pbc(tr_vc::Vector, pbc::Vector) = (tr_vc, pbc)
_tr_and_pbc(tr_vc::Vector, pbc::Bool) = (tr_vc, fill(pbc, length(tr_vc)))
_tr_and_pbc(pbc::Vector) = (zeros(Int, length(pbc)), pbc)
_tr_and_pbc(pbc::Bool) = ([0], [pbc])
function _tr_and_pbc(axis::Int, pbc::Vector)
    tr_vc = zeros(Int, length(pbc))
    tr_vc[axis] = 1
    return (tr_vc, pbc)
end
function _tr_and_pbc(axis::Int, pbc::Bool)
    tr_vc = zeros(Int, axis)
    tr_vc[axis] = 1
    return (tr_vc, fill(pbc, axis))
end
_get_site_indices(::T) where T = error("Cannot convert object of type $T to hopping indices")
_get_site_indices(i::Int) = (i, i)
_get_site_indices(t::NTuple{2, Int}) = t

"""
    hopping([hop_operator,] kwargs...)

A convenient constructor for a `Hopping` object. `hop_operator` can be either a matrix or a number
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
function hopping(hop_operator=1; site_indices=1, pbc=false, kw...)
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
    Hopping(site_indices, tr_vc, pbc, _wrap_operator(hop_operator))
end

==(h1::Hopping, h2::Hopping) =
    all(getproperty(h1, fn) == getproperty(h2, fn) for fn in fieldnames(Hopping))

function show(io::IO, m::MIME"text/plain", hop::Hopping)
    println(io, "Hopping")
    println(io, "Connects site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.translate_uc)")
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
    print(io, "Hopping operator matrix: ")
    show(io, m, hop.hop_operator)
end
dims(h::Hopping) = length(h.translate_uc)
dims_internal(h::Hopping) = size(h.hop_operator)[1]

"""
    promote_dims!(h::Hopping, ndims::Int)

Changes dimension count of `hopping` to `ndims` if possible.
"""
function promote_dims!(h::Hopping, ndims::Int)
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

function _hopping_dest(l::Lattice, hop::Hopping, i_site::LatticeSite)
    i_site.basis_index != hop.site_indices[1] && return nothing
    new_uc = i_site.unit_cell + hop.translate_uc
    for i in 1:dims(l)
        hop.pbc[i] || 1 ≤ new_uc[i] ≤ size(l)[i] || return nothing
    end
    LatticeSite(mod.(new_uc .- 1, l.lattice_size) .+ 1, hop.site_indices[2])
end

check_lattice_fits(::Any, ::Lattice) = nothing
@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
function _hopping_operator!(lop::LatticeOperator, selector, hop::Hopping, field::AbstractField)
    l = lattice(lop)
    check_lattice_fits(selector, l)
    d = dims(l)
    promote_dims!(hop, d)
    trv = SVector{d}(hop.translate_uc)
    for (i, site1) in enumerate(l)
        site2 = _hopping_dest(l, hop, site1)
        site2 === nothing && continue
        j = site_index(l, site2)
        j === nothing && continue
        !_get_bool_value(selector, l, site1, site2) && continue

        p1 = site_coords(l, site1)
        p2 = p1 + trv
        pmod = exp(2π * im * path_integral(field, p1, p2))
        !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
        lop[i, j] = hop.hop_operator * pmod + @view lop[i, j]
        lop[j, i] = (@view lop[i, j])'
    end
    lop
end

@doc raw"""
    hopping_operator([f, ]lattice::Lattice, hopping::Hopping[, field::AbstractField])

Creates a hopping operator:
$$\hat{A} = \sum_{pairs} \hat{t} \hat{c}^\dagger_j \hat{c}_i + h. c.$$

Arguments:
- `f`: takes a `Lattice` and two `LatticeSite`s, returns whether this pair should be included.
- `lattice`: the lattice to create the operator on.
- `hopping`: the `Hopping` object describing the site pairs and the $\hat{t}$ operator.
- `field`: the `AbstractField` object that defines the magnetic field to generate phase factors using Peierls substitution.
"""
function hopping_operator(lf::Function, l::Lattice, hop::Hopping, field::AbstractField=NoField())
    lop = _zero_on_basis(l, hop.hop_operator)
    _hopping_operator!(lop, lf, hop, field)
end
function hopping_operator(l::Lattice, hop::Hopping, field::AbstractField=NoField())
    lop = _zero_on_basis(l, hop.hop_operator)
    _hopping_operator!(lop, nothing, hop, field)
end

"""
    AbstractPairSelector <: Function

A function-like object that accepts a `Lattice` and two `LatticeSite`s and returns whether the site pair is *selected* or not.
Implements a built-in sanity check algorithm to make sure the pair set was defined on a correct lattice.

Define these functions for all subclasses:
- `LatticeModels.lattice(::YourSelector)` must return the lattice your selector was defined on.
- `LatticeModels.match(::YourSelector, site1::LatticeSite, site2::LatticeSite)` must return whether the `(site1, site2)` pair is *selected*.
"""
abstract type AbstractPairSelector <: Function end
match(::T, ::LatticeSite, ::LatticeSite) where {T<:AbstractPairSelector} =
    error("match(::$T, ::LatticeSite, ::LatticeSite) must be explicitly implemented")
lattice(::T) where {T<:AbstractPairSelector} = error("lattice(::$T) must be explicitly implemented")
check_lattice_fits(ps::AbstractPairSelector, l::Lattice) = check_is_sublattice(l, lattice(ps))
(ps::AbstractPairSelector)(::Lattice, site1::LatticeSite, site2::LatticeSite) = match(ps, site1, site2)

"""
    DomainsSelector(domains::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `domains` is the same on two sites.
"""
struct DomainsSelector <: AbstractPairSelector
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
struct PairLhsSelector <: AbstractPairSelector
    lhs::LatticeValue
end
lattice(ps::PairLhsSelector) = lattice(ps.lhs)
match(ps::PairLhsSelector, site1::LatticeSite, ::LatticeSite) =
    ps.lhs[site1]

"""
    PairRhsSelector(rhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `rhs` is true on the first site of the pair.
"""
struct PairRhsSelector <: AbstractPairSelector
    rhs::LatticeValue
end
lattice(ps::PairRhsSelector) = lattice(ps.rhs)
match(ps::PairRhsSelector, ::LatticeSite, site2::LatticeSite) =
    ps.rhs[site2]

macro hopping_operator(for_loop::Expr)
    if for_loop.head !== :for
        throw(ArgumentError("expression must be a for loop, not $(for_loop.head)"))
    elseif length(for_loop.args) != 2
        throw(ArgumentError("malformed for loop")) # This should never happen, but still...
    end
    itr::Expr, body::Expr = for_loop.args
    if !Meta.isexpr(itr, :(=), 2)
        throw(ArgumentError("invalid for loop iteration specification; must be a simple assignment"))
    end
    itr_vars, lattice_var = itr.args
    if !Meta.isexpr(itr_vars, :tuple, 2)
        throw(ArgumentError("invalid for loop iterator variable; must be a 2-tuple"))
    end
    site1_var, site2_var = itr_vars.args
    while body.args[end] isa LineNumberNode
        pop!(body.args)
    end
    dump(body)
    quote
        i = 0
        l = $(esc(lattice_var))
        local matrix = nothing
        local N = 0
        for $(esc(site1_var)) in l
            i += 1
            j = 0
            for $(esc(site2_var)) in l
                j += 1
                if i ≥ j
                    continue
                end
                block_res = $(esc(body))
                if block_res !== nothing
                    if matrix === nothing
                        if size(block_res) == ()
                            N = 1
                        else
                            N = size(block_res)[1]
                        end
                        matrix = zeros(ComplexF64, N * length(l), N * length(l))
                    end
                    matrix[N*(i-1)+1:N*i, N*(j-1)+1:N*j] .= block_res
                    matrix[N*(j-1)+1:N*j, N*(i-1)+1:N*i] .= block_res'
                end
            end
        end
        LatticeArray(Basis(l, N), matrix)
    end
end

"""
    PairSet{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`PairSet`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `PairSet` which connects sites that were connected by `≤n` bonds of the previous `PairSet`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct PairSet{LT<:Lattice} <: AbstractPairSelector
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
function bonds(l::Lattice, hops::Hopping...)
    bs = PairSet(l)
    for h in hops
        promote_dims!(h, dims(l))
    end
    for i in 1:length(l)
        site1 = bs.lattice[i]
        for hop in hops
            site2 = _hopping_dest(l, hop, site1)
            site2 === nothing && continue
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

^(bs1::PairSet, n::Int) =
    PairSet(bs1.lattice, bs1.bmat^n .!= 0)

!(bs::PairSet) =
    PairSet(bs.lattice, Matrix(I, length(bs.lattice), length(bs.lattice)) .| .!(bs.bmat))

function show(io::IO, m::MIME"text/plain", bs::PairSet)
    println(io, "PairSet with $(count(bs.bmat)) bonds")
    print(io, "on ")
    show(io, m, bs.lattice)
end

@recipe function f(bs::PairSet)
    aspect_ratio := :equal
    l = bs.lattice
    pts = Tuple{Float64,Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for i in 1:length(l)
        site1 = bs.lattice[i]
        A = site_coords(l, site1)
        for j in 1:length(l)
            if i != j && bs.bmat[i, j]
                site2 = bs.lattice[j]
                B = site_coords(l, site2)
                T = radius_vector(l, site2, site1)
                push!(pts, Tuple(A))
                push!(pts, Tuple(A + T / 2))
                push!(pts, br_pt)
                push!(pts, Tuple(B))
                push!(pts, Tuple(B - T / 2))
                push!(pts, br_pt)
            end
        end
    end
    label := nothing
    pts
end
