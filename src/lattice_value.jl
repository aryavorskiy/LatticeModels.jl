using LinearAlgebra, Statistics, Logging
import Base: length, size, getindex, setindex!, eltype, copyto!, show, ==

struct LatticeValue{T}
    lattice::AbstractLattice
    vector::Vector{T}
    LatticeValue(lattice::AbstractLattice, vector::Vector{T}) where {T} =
        new{T}(lattice, vector)
end

LatticeValue(f, l::AbstractLattice) =
    (lf = _propagate_lattice_args(f, l); LatticeValue(l, [lf(l, site) for site in l]))

==(lv1::LatticeValue, lv2::LatticeValue) = (lv1.lattice == lv2.lattice) && (lv1.vector == lv2.vector)
eltype(::LatticeValue{T}) where {T} = T
length(lv::LatticeValue) = length(lv.vector)
size(lv::LatticeValue) = size(lv.vector)
iterate(lv::LatticeValue) = iterate(lv.vector)
iterate(lv::LatticeValue, state) = iterate(lv.vector, state)
getindex(lv::LatticeValue, cartesian_i::CartesianIndex{1}) = lv.vector[cartesian_i]

struct LVStyle <: Broadcast.BroadcastStyle end
copyto!(lv::LatticeValue, src::Broadcast.Broadcasted{LVStyle}) = (copyto!(lv.vector, src); return lv)
Base.broadcastable(lv::LatticeValue) = lv
Base.BroadcastStyle(::Type{<:LatticeValue}) = LVStyle()
Base.BroadcastStyle(bs::Broadcast.BroadcastStyle, ::LVStyle) =
    error("cannot broadcast LatticeValue along style $bs")
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LVStyle) = LVStyle()

function Base.similar(bc::Broadcast.Broadcasted{LVStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc.args)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(args::Tuple) = _extract_lattice(args[begin], Base.tail(args))
_extract_lattice(::Tuple{}) = nothing
_extract_lattice(::Any, rem_args::Tuple) = _extract_lattice(rem_args)
_extract_lattice(lv::LatticeValue, rem_args::Tuple) = _extract_lattice_s(lv.lattice, rem_args)

_extract_lattice_s(l::AbstractLattice, args::Tuple) =
    _extract_lattice_s(l, args[begin], Base.tail(args))
_extract_lattice_s(l::AbstractLattice, ::Tuple{}) = l
_extract_lattice_s(l::AbstractLattice, ::Any, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)
function _extract_lattice_s(l::AbstractLattice, lv::LatticeValue, rem_args::Tuple)
    lv.lattice != l && throw(ArgumentError("lattice mismatch:\n$l\n$(lv.lattice)"))
    _extract_lattice_s(l, rem_args)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeValue{T}) where {T}
    println(io, "LatticeValue with inner type $T")
    print(io, "on ")
    show(io, m, lv.lattice)
end

_supports_heatmap(::FiniteBravaisLattice) = false
_supports_heatmap(::SquareLattice) = true
_supports_heatmap(sl::SubLattice) = _supports_heatmap(sl.lattice)
_heatmap_axes(l::SquareLattice) = [-(ax - 1)/2:(ax-1)/2 for ax in _sz(l)]
_heatmap_axes(sl::SubLattice) = _heatmap_axes(sl.lattice)
_heatmap_vals(::FiniteBravaisLattice, vec) = vec
function _heatmap_vals(sl::SubLattice, vec)
    @assert length(sl) == length(vec)
    i = 1
    len = length(sl.mask)
    newvals = fill(NaN, len)
    @inbounds for j in 1:len
        if sl.mask[j]
            newvals[j] = vec[i]
            i += 1
        end
    end
    newvals
end

@recipe function f(lv::LatticeValue)
    _supports_heatmap(lv.lattice) && (seriestype --> :heatmap)
    if _supports_heatmap(lv.lattice) && plotattributes[:seriestype] == :heatmap
        b = _heatmap_vals(lv.lattice, lv.vector)
        aspect_ratio := :equal
        _heatmap_axes(lv.lattice)..., reshape(b, Tuple(_sz(lv.lattice)))
    else
        lv.lattice, lv.vector
    end
end

# TODO: add logical indexing support
