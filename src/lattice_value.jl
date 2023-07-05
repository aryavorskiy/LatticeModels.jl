using LinearAlgebra, Logging
import Base: ==

struct LatticeValueWrapper{VT<:AbstractVector,LatticeSym}
    lattice::Lattice{LatticeSym}
    values::VT
    function LatticeValueWrapper(lattice::Lattice{LatticeSym}, values::VT) where {VT,LatticeSym}
        length(lattice) != length(values) &&
            throw(DimensionMismatch("vector has length $(length(values)), lattice has length $(length(lattice))"))
        new{VT,LatticeSym}(lattice, values)
    end
end

lattice(lvw::LatticeValueWrapper) = lvw.lattice

Base.copy(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), copy(lvw.values))
Base.length(lvw::LatticeValueWrapper) = length(lvw.values)
Base.size(lvw::LatticeValueWrapper) = size(lvw.values)
function _to_index(lvw::LatticeValueWrapper, site::LatticeSite)
    i = CartesianIndex(site_index(lattice(lvw), site))
    i === nothing ? throw(BoundsError(lvw, site)) : return i
end
_to_index(::LatticeValueWrapper, i::CartesianIndex{1}) = i
Base.getindex(lvw::LatticeValueWrapper, i) = getindex(lvw.values, _to_index(lvw, i))
Base.setindex!(lvw::LatticeValueWrapper, val, i) = setindex!(lvw.values, val, _to_index(lvw, i))
Base.eltype(lvw::LatticeValueWrapper) = eltype(lvw.values)
Base.eachindex(lvw::LatticeValueWrapper) = lattice(lvw)
Base.iterate(lvw::LatticeValueWrapper, s...) = iterate(lvw.values, s...)
Base.pairs(lvw::LatticeValueWrapper) = Iterators.map(=>, lvw.lattice, lvw.values)

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
LatticeValue(lf, l::Lattice) = LatticeValue(l, [lf(site) for site in l])

"""
    coord_values(l::Lattice)

Generates a tuple of `LatticeValue`s representing spatial coordinates.
"""
coord_values(l::Lattice) =
    [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]

Base.rand(l::Lattice) = LatticeValue(l, rand(length(l)))
Base.rand(T::Type, l::Lattice) = LatticeValue(l, rand(T, length(l)))
Base.randn(l::Lattice) = LatticeValue(l, randn(length(l)))
Base.randn(T::Type, l::Lattice) = LatticeValue(l, randn(T, length(l)))
Base.fill(value, l::Lattice) = LatticeValue(l, fill(value, length(l)))
Base.fill!(lv::LatticeValue, value) = (fill!(lv.values, value); lv)
Base.zero(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), zero(lvw.values))
Base.zeros(l::Lattice) = fill(0., l)
Base.zeros(T::Type, l::Lattice) = fill(zero(T), l)
Base.one(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), one(lvw.values))
Base.ones(l::Lattice) = fill(1., l)
Base.ones(T::Type, l::Lattice) = fill(one(T), l)

==(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper) = (lvw1.lattice == lvw2.lattice) && (lvw1.values == lvw2.values)

struct LVWStyle <: Broadcast.BroadcastStyle end
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{LVWStyle}) = (copyto!(lvw.values, src); return lvw)
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{Broadcast.DefaultArrayStyle{0}}) = (copyto!(lvw.values, src); return lvw)
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

function Base.show(io::IO, mime::MIME"text/plain", lv::LatticeValueWrapper)
    print(io, "$(typeof(lv)) on a ")
    show(io, mime, lattice(lv))
    if !get(io, :compact, false)
        print(io, "\nValues stored in a $(typeof(lv.values)):\n")
        Base.show_vector(IOContext(io, :compact => true), lv.values)
    end
end

Base.@propagate_inbounds function Base.getindex(l::Lattice{LatticeSym,N,NB},
        lv_mask::LatticeValue{Bool,LatticeSym}) where {LatticeSym,N,NB}
    @boundscheck check_is_sublattice(l, lattice(lv_mask))
    new_mask = zero(l.mask)
    new_mask[lv_mask.lattice.mask] = lv_mask.values
    Lattice(LatticeSym, size(l), bravais(l), vec(new_mask .& l.mask))
end

Base.@propagate_inbounds function Base.getindex(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    new_l  = lattice(lv)[lv_mask]
    LatticeValueWrapper(new_l, lv.values[new_l.mask[lv.lattice.mask]])
end

Base.@propagate_inbounds function Base.Broadcast.dotview(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    new_l  = lattice(lv)[lv_mask]
    LatticeValueWrapper(new_l, view(lv.values, new_l.mask[lv.lattice.mask]))
end
Base.Broadcast.dotview(lv::LatticeValueWrapper; kw...) = Base.Broadcast.dotview(lv, _kws_to_mask(lattice(lv), kw))

Base.@propagate_inbounds function Base.setindex!(lv::LatticeValueWrapper, lv_rhs::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck begin
        check_is_sublattice(lattice(lv), lattice(lv_mask))
        check_is_sublattice(lattice(lv), lattice(lv_rhs))
    end
    new_mask = zero(lv.lattice.mask)
    new_mask[lv_rhs.lattice.mask] = lv_mask.values
    lv.values[new_mask[lv.lattice.mask]] = lv_rhs.values[new_mask[lv_rhs.lattice.mask]]
end

function _parse_axis_descriptor(coord::Symbol)
    coord === :x && return (:c, 1)
    coord === :y && return (:c, 2)
    coord === :z && return (:c, 3)
    axis_s = string(coord)
    axis = tryparse(Int, axis_s[2:end])
    if axis isa Int && axis > 0
        axis_s[1] == 'x' && return (:c, axis)
        axis_s[1] == 'j' && return (:j, axis)
    end
    error("failed to parse axis descriptor '$coord'")
end
function _parse_axis_descriptor(i::Int)
    i ≤ 0 && error("axis number must be a positive integer, not $i")
    (:c, i)
end
function _kws_to_mask(l, @nospecialize(kw))
    coords = collect_coords(l)
    relmask = fill(true, length(l))
    for (descr, val) in kw
        if descr === :index
            @. relmask &= (site.basis_index == val for site in l)
            continue
        end

        type, axis_no = _parse_axis_descriptor(descr)
        axis_no > dims(l) && error(
            "axis number $axis_no (descriptor '$descr') exceeds lattice dimensionality $(dims(l))")
        type === :j && @.(relmask &= (site.unit_cell[axis_no] == val for site in l))
        type === :c && @.(relmask &= abs(coords[axis_no, :] - val) < eps())
    end
    cnt = count(relmask)
    cnt == 1 ? CartesianIndex(findfirst(relmask)) : LatticeValue(l, relmask)
end
Base.getindex(lvw::LatticeValueWrapper; kw...) = getindex(lvw, _kws_to_mask(lattice(lvw), kw))
Base.getindex(l::Lattice; kw...) = getindex(l, _kws_to_mask(l, kw))
Base.setindex!(lvw::LatticeValueWrapper, rhs; kw...) = setindex!x(lvw, rhs, _kws_to_mask(lattice(lvw), kw))

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

"""
    project(lv::LatticeValue, axis)

Creates a mapping from site coordinates to values of `lv`.
The coordinate axis to project the sites onto can be set with the `axis` argument -
it can be either an integer from 1 to 3 or an axis descriptor `Symbol`.
"""
function project(lv::PlottableLatticeValue, axis)
    type, axis_no = _parse_axis_descriptor(axis)
    l = lattice(lv)
    axis_no > dims(lattice(lv)) && error(
        "axis number $axis_no (descriptor '$axis') exceeds lattice dimensionality $(dims(l))")
    crds = collect_coords(l)
    pr_crds = type === :c ? pr_crds = crds[axis_no, :] :
        crds' * normalize(bravais(l).translation_vectors[:, axis_no])
    perm = sortperm(pr_crds)
    pr_crds[perm], lv.values[perm]
end
