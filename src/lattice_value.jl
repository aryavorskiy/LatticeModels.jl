using LinearAlgebra, Statistics, Logging
import Base: length, size, getindex, setindex!, eltype, copyto!, show, ==

"""
    LatticeValue{T, LatticeSym}

Represents a value of type `T` on a `Lattice{LatticeSym}` lattice.

Fields:
- lattice: the `Lattice` object the value is defined on
- values: the values on different sites
"""
struct LatticeValue{T,LatticeSym}
    lattice::Lattice{LatticeSym}
    values::Vector{T}
    """
        LatticeValue(lattice::Lattice, vector::AbstractVector)

    Constructs a LatticeValue object.
    """
    function LatticeValue(lattice::Lattice{LatticeSym}, values::AbstractVector{T}) where {T,LatticeSym}
        length(lattice) != length(values) &&
            error("inconsistent vector length:\nlattice: $(length(lattice)), vector: $(length(values))")
        new{T,LatticeSym}(lattice, values)
    end
end

LatticeValue(lf, l::Lattice) = LatticeValue(l, [lf(site, site_coords(l, site)) for site in l])

lattice(l::LatticeValue) = l.lattice

"""
    coord_values(lattice::Lattice)

Generates a tuple of `LatticeValue`s representing coordinate functions.
"""
coord_values(l::Lattice) = [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]

import Base: rand, randn
rand(l::Lattice) = LatticeValue(l, rand(length(l)))
rand(T::Type, l::Lattice) = LatticeValue(l, rand(T, length(l)))
randn(l::Lattice) = LatticeValue(l, randn(length(l)))
randn(T::Type, l::Lattice) = LatticeValue(l, randn(T, length(l)))

==(lv1::LatticeValue, lv2::LatticeValue) = (lv1.lattice == lv2.lattice) && (lv1.values == lv2.values)
eltype(::LatticeValue{T}) where {T} = T
length(lv::LatticeValue) = length(lv.values)
size(lv::LatticeValue) = size(lv.values)
getindex(lv::LatticeValue, cartesian_i::CartesianIndex{1}) = lv.values[cartesian_i]

struct LVStyle <: Broadcast.BroadcastStyle end
copyto!(lv::LatticeValue, src::Broadcast.Broadcasted{LVStyle}) = (copyto!(lv.values, src); return lv)
Base.broadcastable(lv::LatticeValue) = lv
Base.BroadcastStyle(::Type{<:LatticeValue}) = LVStyle()
Base.BroadcastStyle(bs::Broadcast.BroadcastStyle, ::LVStyle) =
    error("cannot broadcast LatticeValue along style $bs")
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LVStyle) = LVStyle()

function Base.similar(bc::Broadcast.Broadcasted{LVStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(bc::Broadcast.Broadcasted) = _extract_lattice(bc.args)
_extract_lattice(lv::LatticeValue) = lv.lattice
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
    l2 != l && throw(ArgumentError("lattice mismatch:\n$l\n$l2"))
    _extract_lattice_s(l, rem_args)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeValue{T}) where {T}
    println(io, "LatticeValue with inner type $T")
    print(io, "on ")
    show(io, m, lattice(lv))
end

function getindex(l::Lattice{LatticeSym,N,NB}, lvm::LatticeValue{Bool,LatticeSym}) where {LatticeSym,N,NB}
    l != lattice(lvm) && error("lattice mismatch")
    new_mask = zero(l.mask)
    new_mask[l.mask] = lvm.values
    Lattice(LatticeSym, size(l), bravais(l), new_mask)
end

function getindex(lv::LatticeValue, lvm::LatticeValue{Bool})
    lattice(lv) != lattice(lvm) && error("lattice mismatch")
    LatticeValue(lv.lattice[lvm], lv.values[lvm.values])
end

Base.@propagate_inbounds function getindex(lv::LatticeValue, site::LatticeSite)
    i = site_index(site, lattice(lv))
    @boundscheck i === nothing && throw(BoundsError(lv, site))
    lv.values[i]
end

_heatmap_axes(l::SquareLattice) = [-(ax - 1)/2:(ax-1)/2 for ax in size(l)]

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
    project(lattice_value, axis)

Creates a mapping from site coordinates to values of `lattice_value`.
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
