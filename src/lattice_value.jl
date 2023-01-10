using LinearAlgebra, Statistics, Logging
import Base: length, size, pairs, getindex, setindex!, eachindex, eltype, copyto!, show, ==

struct LatticeValueWrapper{VT<:AbstractVector,LatticeSym}
    lattice::Lattice{LatticeSym}
    values::VT
    function LatticeValueWrapper(lattice::Lattice{LatticeSym}, values::VT) where {VT,LatticeSym}
        length(lattice) != length(values) &&
            error("inconsistent vector length:\nlattice: $(length(lattice)), vector: $(length(values))")
        new{VT,LatticeSym}(lattice, values)
    end
end

lattice(lvw::LatticeValueWrapper) = lvw.lattice
lattice(l::Lattice) = l

size(lvw::LatticeValueWrapper) = size(lvw.values)
length(lvw::LatticeValueWrapper) = length(lvw.values)
function _to_index(lvw::LatticeValueWrapper, site::LatticeSite)
    i = CartesianIndex(site_index(lattice(lvw), site))
    i === nothing ? error("index conversion failed") : return i
end
_to_index(::LatticeValueWrapper, i::CartesianIndex{1}) = i
getindex(lvw::LatticeValueWrapper, i) = getindex(lvw.values, _to_index(lvw, i))
setindex!(lvw::LatticeValueWrapper, val, i) =
    setindex!(lvw.values, val, _to_index(lvw, i))
eltype(lvw::LatticeValueWrapper) = eltype(lvw.values)
eachindex(lvw::LatticeValueWrapper) = lattice(lvw)
iterate(lvw::LatticeValueWrapper, s...) = iterate(lvw.values, s...)
pairs(lvw::LatticeValueWrapper) = Iterators.map(=>, lvw.lattice, lvw.values)

"""
    LatticeValue{T, LatticeSym}

Represents a value of type `T` on a `Lattice{LatticeSym}` lattice.

Fields:
- lattice: the `Lattice` object the value is defined on
- values: the values on different sites
"""
const LatticeValue{T, LT} = LatticeValueWrapper{Vector{T}, LT}

"""
    LatticeValue(l::Lattice, v::AbstractVector)

Constructs a LatticeValue object.
"""
LatticeValue(l::Lattice, v::AbstractVector) = LatticeValueWrapper(l, convert(Vector, v))
LatticeValue(lf, l::Lattice) = LatticeValue(l, [lf(site, site_coords(l, site)) for site in l])

"""
    coord_values(l::Lattice)

Generates a tuple of `LatticeValue`s representing spatial coordinates.
"""
coord_values(l::Lattice) = [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]

import Base: rand, randn, fill, fill!, zero, zeros, one, ones
rand(l::Lattice) = LatticeValue(l, rand(length(l)))
rand(T::Type, l::Lattice) = LatticeValue(l, rand(T, length(l)))
randn(l::Lattice) = LatticeValue(l, randn(length(l)))
randn(T::Type, l::Lattice) = LatticeValue(l, randn(T, length(l)))
fill(value, l::Lattice) = LatticeValue(l, fill(value, length(l)))
fill!(lv::LatticeValue, value) = (fill!(lv.values, value); lv)
zero(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), zero(lvw.values))
zeros(l::Lattice) = fill(0., l)
zeros(T::Type, l::Lattice) = fill(zero(T), l)
one(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), one(lvw.values))
ones(l::Lattice) = fill(1., l)
ones(T::Type, l::Lattice) = fill(one(T), l)

==(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper) = (lvw1.lattice == lvw2.lattice) && (lvw1.values == lvw2.values)

struct LVWStyle <: Broadcast.BroadcastStyle end
copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{LVWStyle}) = (copyto!(lvw.values, src); return lvw)
copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{Broadcast.DefaultArrayStyle{0}}) = (copyto!(lvw.values, src); return lvw)
Base.broadcastable(lvw::LatticeValueWrapper) = lvw
Base.BroadcastStyle(::Type{<:LatticeValueWrapper}) = LVWStyle()
Base.BroadcastStyle(bs::Broadcast.BroadcastStyle, ::LVWStyle) =
    error("cannot broadcast LatticeValue along style $bs")
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LVWStyle) = LVWStyle()

function Base.similar(bc::Broadcast.Broadcasted{LVWStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(bc::Broadcast.Broadcasted) = _extract_lattice(bc.args)
_extract_lattice(lv::LatticeValueWrapper) = lv.lattice
_extract_lattice(x) = x
_extract_lattice(::Tuple{}) = nothing
_extract_lattice(args::Tuple) =
    _extract_lattice(_extract_lattice(args[begin]), Base.tail(args))
_extract_lattice(::Any, rem_args::Tuple) = _extract_lattice(rem_args)
_extract_lattice(l::Lattice, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)

_extract_lattice_s(l::Lattice, args::Tuple) =
    _extract_lattice_s(l, _extract_lattice(args[begin]), Base.tail(args))
_extract_lattice_s(l::Lattice, ::Tuple{}) = l
_extract_lattice_s(l::Lattice, ::Any, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)
function _extract_lattice_s(l::Lattice, l2::Lattice, rem_args::Tuple)
    check_lattice_match(l, l2)
    _extract_lattice_s(l, rem_args)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeValue{T}) where {T}
    println(io, "LatticeValue with eltype $T\non ")
    show(io, m, lattice(lv))
end
function show(io::IO, m::MIME"text/plain", lv::LatticeValueWrapper{VT}) where {VT}
    println(io, "LatticeValueWrapper with inner type $VT\non ")
    show(io, m, lattice(lv))
end

Base.@propagate_inbounds function getindex(l::Lattice{LatticeSym,N,NB},
        lv_mask::LatticeValue{Bool,LatticeSym}) where {LatticeSym,N,NB}
    @boundscheck check_is_sublattice(l, lattice(lv_mask))
    new_mask = zero(l.mask)
    new_mask[lv_mask.lattice.mask] = lv_mask.values
    Lattice(LatticeSym, size(l), bravais(l), vec(new_mask .& l.mask))
end

Base.@propagate_inbounds function getindex(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    new_l  = lattice(lv)[lv_mask]
    LatticeValueWrapper(new_l, lv.values[new_l.mask[lv.lattice.mask]])
end

Base.@propagate_inbounds function Base.maybeview(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    new_l  = lattice(lv)[lv_mask]
    LatticeValueWrapper(new_l, view(lv.values, new_l.mask[lv.lattice.mask]))
end

Base.@propagate_inbounds function setindex!(lv::LatticeValueWrapper, lv_rhs::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck begin
        check_is_sublattice(lattice(lv), lattice(lv_mask))
        check_is_sublattice(lattice(lv), lattice(lv_rhs))
    end
    new_mask = zero(lv.lattice.mask)
    new_mask[lv_rhs.lattice.mask] = lv_mask.values
    lv.values[new_mask[lv.lattice.mask]] = lv_rhs.values[new_mask[lv_rhs.lattice.mask]]
end

raw"""
    macro_cell_values(lv::LatticeValue)

Returng an array of the values of `lv` on its macro cell.
The $i$-th element of the array corresponds to the $i$-th site of the macro cell.
If the element is `NaN`, it means that the corresponding site is not present in the `lv`'s lattice.

This function might be quite useful in custom plot recipes.
"""
function macro_cell_values(lv::LatticeValue{<:Number})
    i = 1
    len = length(lv.lattice.mask)
    newvals = fill(NaN, len)
    @inbounds for j in 1:len
        if lv.lattice.mask[j]
            newvals[j] = lv.values[i]
            i += 1
        end
    end
    newvals
end

"""
    plot_fallback(lv::LatticeValue)

Creates a copy of `lv` lattice value with its `LatticeSym` overwritten to `:plot_fallback`.
Use it to invoke the default plot recipe for `LatticeValues` when defining a custom one.
"""
function plot_fallback(lv::LatticeValue)
    l = lattice(lv)
    new_l = Lattice(:plot_fallback, size(l), bravais(l), l.mask)
    LatticeValue(new_l, lv.values)
end

const PlottableLatticeValue{LT} = LatticeValue{<:Number, LT}

@recipe function f(lv::PlottableLatticeValue{:square})
    seriestype --> :heatmap
    if plotattributes[:seriestype] === :heatmap
        aspect_ratio := :equal
        axes_lims = [-(ax - 1)/2:(ax-1)/2 for ax in size(lattice(lv))]
        heatmap_values = reshape(macro_cell_values(lv), reverse(size(lattice(lv))))'
        axes_lims..., heatmap_values
    else
        plot_fallback(lv)
    end
end

@recipe function f(lv::PlottableLatticeValue)
    lv.lattice, lv.values
end

"""
    project(lv::LatticeValue, axis)

Creates a mapping from site coordinates to values of `lv`.
The coordinate axis to project the sites onto can be set with the `axis` argument -
it can be either an integer from 1 to 3 or a `Symbol` (`:x`, `:y` or `:z`).
"""
function project(lv::PlottableLatticeValue, axis)
    axis_no = axis isa Int ? axis :
        axis === :x ? 1 :
        axis === :y ? 2 :
        axis === :z ? 3 : 0
    axis_no ∉ 1:3 && error("unsupported projection axis '$axis'")
    crds = collect_coords(lattice(lv))[axis_no, :]
    perm = sortperm(crds)
    crds[perm], lv.values[perm]
end
