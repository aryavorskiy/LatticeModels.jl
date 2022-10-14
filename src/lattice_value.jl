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

LatticeValue(lf, l::Lattice) = LatticeValue(l, [lf(site, coords(l, site)) for site in l])

"""
    coord_values(lattice::Lattice)

Generates a tuple of `LatticeValue`s representing coordinate functions.
"""
coord_values(l::Lattice) = [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]

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
    show(io, m, lv.lattice)
end

function getindex(l::Lattice{LatticeSym,N,NB}, lvm::LatticeValue{Bool,LatticeSym}) where {LatticeSym,N,NB}
    l != lvm.lattice && error("lattice mismatch")
    new_mask = zero(l.mask)
    new_mask[l.mask] = lvm.values
    Lattice(LatticeSym, size(l), bravais(l), new_mask)
end

function getindex(lv::LatticeValue, lvm::LatticeValue{Bool})
    lv.lattice != lvm.lattice && error("lattice mismatch")
    LatticeValue(lv.lattice[lvm], lv.values[lvm.values])
end

_heatmap_axes(l::SquareLattice) = [-(ax - 1)/2:(ax-1)/2 for ax in size(l)]
function _heatmap_vals(slv::LatticeValue{<:Number,:square})
    i = 1
    len = length(slv.lattice.mask)
    newvals = fill(NaN, len)
    @inbounds for j in 1:len
        if slv.lattice.mask[j]
            newvals[j] = slv.values[i]
            i += 1
        end
    end
    newvals
end

"""
    plot_fallback(lv::LatticeValue)

Creates a copy of `lv` lattice value with its `LatticeSym` overwritten to `:uncertain`.
Use it to invoke the default plot recipe for `LatticeValues` when defining a custom one.
"""
function plot_fallback(lv::LatticeValue)
    l = lv.lattice
    new_l = Lattice(:uncertain, size(l), bravais(l), l.mask)
    LatticeValue(new_l, lv.values)
end

@recipe function f(lv::LatticeValue{<:Number,:square})
    if get(plotattributes, :seriestype, :unknown) === :heatmap
        aspect_ratio := :equal
        _heatmap_axes(lv.lattice)..., reshape(_heatmap_vals(lv), size(lv.lattice))
    else
        plot_fallback(lv)
    end
end

@recipe function f(lv::LatticeValue; project_axis=:none)
    if project_axis !== :none
        axis_no =   project_axis isa Int ? project_axis :
                    project_axis === :x ? 1 :
                    project_axis === :y ? 2 :
                    project_axis === :z ? 3 : 0
        axis_no ∉ 1:3 && error("unsupported projection axis '$project_axis'")
        crds = collect_coords(lv.lattice)[axis_no, :]
        perm = sortperm(crds)
        crds[perm], lv.values[perm]
    else
        lv.lattice, lv.values
    end
end
