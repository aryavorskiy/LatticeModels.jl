using LinearAlgebra, Statistics, Logging
import Base: length, size, getindex, setindex!, eltype, copyto!, show, ==

struct LatticeValue{T,LatticeSym}
    lattice::Lattice
    vector::Vector{T}
    function LatticeValue(lattice::Lattice{LatticeSym}, vector::AbstractVector{T}) where {T,LatticeSym}
        length(lattice) != length(vector) &&
            error("inconsistent vector length")
        new{T,LatticeSym}(lattice, vector)
    end
end

LatticeValue(lf, l::Lattice) = LatticeValue(l, [lf(site, coords(l, site)) for site in l])

function coord_values(l::Lattice)
    d = dims(l)
    i = 1
    xyz_values = zeros(length(l), d)
    for site in l
        crd = coords(l, site)
        xyz_values[i, :] = crd
        i += 1
    end
    [LatticeValue(l, vec) for vec in eachcol(xyz_values)]
end

==(lv1::LatticeValue, lv2::LatticeValue) = (lv1.lattice == lv2.lattice) && (lv1.vector == lv2.vector)
eltype(::LatticeValue{T}) where {T} = T
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
    new_mask[l.mask] = lvm.vector
    Lattice(LatticeSym, size(l), bravais(l), new_mask)
end

function getindex(lv::LatticeValue, lvm::LatticeValue{Bool})
    lv.lattice != lvm.lattice && error("lattice mismatch")
    LatticeValue(lv.lattice[lvm], lv.vector[lvm.vector])
end

_heatmap_axes(l::SquareLattice) = [-(ax - 1)/2:(ax-1)/2 for ax in size(l)]
function _heatmap_vals(slv::LatticeValue{<:Number,:square})
    i = 1
    len = length(slv.lattice.mask)
    newvals = fill(NaN, len)
    @inbounds for j in 1:len
        if slv.lattice.mask[j]
            newvals[j] = slv.vector[i]
            i += 1
        end
    end
    newvals
end

@recipe function f(lv::LatticeValue{<:Number,:square})
    seriestype --> :heatmap
    if plotattributes[:seriestype] == :heatmap
        aspect_ratio := :equal
        _heatmap_axes(lv.lattice)..., reshape(_heatmap_vals(lv), size(lv.lattice))
    else
        lv.lattice, lv.vector
    end
end

@recipe function f(lv::LatticeValue)
    lv.lattice, lv.vector
end
