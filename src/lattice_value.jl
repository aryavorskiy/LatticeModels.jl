using LinearAlgebra, Statistics, Logging
import Base: length, size, getindex, setindex!, eltype, copyto!, show, ==

struct LatticeValue{LatticeType, T}
    lattice::Lattice
    vector::Vector{T}
    function LatticeValue(lattice::Lattice{LatticeType}, vector::Vector{T}) where {T, LatticeType}
        length(lattice) != length(vector) &&
            error("inconsistent vector length")
        new{LatticeType, T}(lattice, vector)
    end
end

LatticeValue(f, l::Lattice) =
    (lf = _propagate_lattice_args(f, l); LatticeValue(l, [lf(l, site) for site in l]))

function coord_values(@nospecialize l::Lattice)
    d = dims(l)
    i = 1
    xyz_values = zeros(length(l), d)
    for site in l
        crd = coords(l, site)
        xyz_values[i, :] = crd
        i += 1
    end
    [LatticeValue(l, vec) for vec in eachrow(xyz_values)]
end

function site_indices_values(l::Lattice)
    d = dims(l)
    i = 1
    ijk_values = zeros(length(l), d+1)
    for site::LatticeIndex in l
        crd = coords(l, site)
        xyz_values[i, 1:d] = site.unit_cell
        xyz_values[i, d+1] = site.basis_index
        i += 1
    end
    [LatticeValue(l, vec) for vec in eachrow(ijk_values)]
end

==(lv1::LatticeValue, lv2::LatticeValue) = (lv1.lattice == lv2.lattice) && (lv1.vector == lv2.vector)
eltype(::LatticeValue{LatticeType, T}) where {LatticeType, T} = T
length(lv::LatticeValue) = length(lv.vector)
size(lv::LatticeValue) = size(lv.vector)
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

_extract_lattice_s(l::Lattice, args::Tuple) =
    _extract_lattice_s(l, args[begin], Base.tail(args))
_extract_lattice_s(l::Lattice, ::Tuple{}) = l
_extract_lattice_s(l::Lattice, ::Any, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)
function _extract_lattice_s(l::Lattice, lv::LatticeValue, rem_args::Tuple)
    lv.lattice != l && throw(ArgumentError("lattice mismatch:\n$l\n$(lv.lattice)"))
    _extract_lattice_s(l, rem_args)
end

function show(io::IO, m::MIME"text/plain", lv::LatticeValue{T}) where {T}
    println(io, "LatticeValue with inner type $T")
    print(io, "on ")
    show(io, m, lv.lattice)
end

_heatmap_axes(l::Lattice{:square}) = [-(ax - 1)/2:(ax-1)/2 for ax in size(l)]
function _heatmap_vals(sl::Lattice{:square}, vec)
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

@recipe function f(lv::LatticeValue{:square})
    seriestype --> :heatmap
    if plotattributes[:seriestype] == :heatmap
        b = _heatmap_vals(lv.lattice, lv.vector)
        aspect_ratio := :equal
        _heatmap_axes(lv.lattice)..., reshape(b, size(lv.lattice))
    else
        lv.lattice, lv.vector
    end
end

@recipe function f(lv::LatticeValue)
    lv.lattice, lv.vector
end

# TODO: add logical indexing support
