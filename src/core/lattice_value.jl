using LinearAlgebra, Logging, QuantumOpticsBase

struct LatticeValueWrapper{LT<:Lattice, VT<:AbstractVector}
    lattice::LT
    values::VT
    function LatticeValueWrapper(lattice::LT, values::VT) where {LT,VT}
        length(lattice) != length(values) &&
            throw(DimensionMismatch("vector has length $(length(values)), lattice has length $(length(lattice))"))
        new{LT,VT}(lattice, values)
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
_to_index(::LatticeValueWrapper, i) = i
Base.getindex(lvw::LatticeValueWrapper, i) = getindex(lvw.values, _to_index(lvw, i))
Base.setindex!(lvw::LatticeValueWrapper, val, i) = setindex!(lvw.values, val, _to_index(lvw, i))
Base.eltype(lvw::LatticeValueWrapper) = eltype(lvw.values)
Base.eachindex(lvw::LatticeValueWrapper) = lattice(lvw)
Base.iterate(lvw::LatticeValueWrapper, s...) = iterate(lvw.values, s...)
Base.pairs(lvw::LatticeValueWrapper) = Iterators.map(=>, lvw.lattice, lvw.values)

"""
    LatticeValue{T, LT}

Represents a value of type `T` on a `Lattice{LatticeSym}` lattice.

Fields:
- lattice: the `Lattice` object the value is defined on
- values: the values on different sites
"""
const LatticeValue{T, LT} = LatticeValueWrapper{LT, Vector{T}}

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
coord_value(l::Lattice, a) = LatticeValue(l, [get_coord(site, parse_axis_sym(a)) for site in l])

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

Base.:(==)(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper) = (lvw1.lattice == lvw2.lattice) && (lvw1.values == lvw2.values)

struct LatticeStyle <: Broadcast.BroadcastStyle end
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{LatticeStyle}) = (copyto!(lvw.values, src); return lvw)
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{Broadcast.DefaultArrayStyle{0}}) = (copyto!(lvw.values, src); return lvw)
Base.broadcastable(lvw::LatticeValueWrapper) = lvw
Base.broadcastable(l::Lattice) = l
Base.BroadcastStyle(::Type{<:LatticeValueWrapper}) = LatticeStyle()
Base.BroadcastStyle(::Type{<:Lattice}) = LatticeStyle()
Base.BroadcastStyle(bs::Broadcast.BroadcastStyle, ::LatticeStyle) =
    error("cannot broadcast LatticeValue along style $bs")
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LatticeStyle) = LatticeStyle()

function Base.similar(bc::Broadcast.Broadcasted{LatticeStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(bc::Broadcast.Broadcasted) = _extract_lattice(bc.args)
_extract_lattice(lv::LatticeValueWrapper) = lv.lattice
_extract_lattice(ext::Broadcast.Extruded) = _extract_lattice(ext.x)
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
    check_samelattice(l, l2)
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

pointer_inds(l::LT, lv_mask::LatticeValue{Bool, LT}) where LT<:Lattice =
    Int[site_index(l, lattice(lv_mask).pointers[i])
        for i in eachindex(lv_mask.values) if lv_mask.values[i]]
Base.@propagate_inbounds function Base.getindex(l::Lattice, lv_mask::LatticeValue{Bool})
    @boundscheck check_issublattice(lattice(lv_mask), l)
    inds = pointer_inds(l, lv_mask)
    Lattice(l.bravais, l.pointers[inds])
end

Base.@propagate_inbounds function Base.getindex(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck check_samelattice(lv, lv_mask)
    inds = pointer_inds(lattice(lv), lv_mask)
    LatticeValueWrapper(lattice(lv)[inds], lv.values[inds])
end

Base.@propagate_inbounds function Base.Broadcast.dotview(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck check_samelattice(lv, lv_mask)
    inds = pointer_inds(lattice(lv), lv_mask)
    LatticeValueWrapper(lattice(lv)[inds], @view lv.values[inds])
end
Base.Broadcast.dotview(lv::LatticeValueWrapper; kw...) = Base.Broadcast.dotview(lv, _kws_to_mask(lattice(lv), kw))

Base.@propagate_inbounds function Base.setindex!(lv::LatticeValueWrapper, lv_rhs::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck begin
        check_issublattice(lattice(lv_mask), lattice(lv))
        check_issublattice(lattice(lv_rhs), lattice(lv))
    end
    inds_l = pointer_inds(lattice(lv), lv_mask)
    inds_r = pointer_inds(lattice(lv_rhs), lv_mask)
    lv.values[inds_l] = lv_rhs.values[inds_r]
end

function _kws_to_mask(l, @nospecialize(kw))
    relmask = fill(true, length(l))
    for (descr, val) in kw
        axis = parse_axis_sym(descr)
        relmask .&= (get_coord(site, axis) in val for site in l)
    end
    cnt = count(relmask)
    cnt == 1 ? CartesianIndex(findfirst(relmask)) : LatticeValue(l, relmask)
end
Base.getindex(lvw::LatticeValueWrapper; kw...) = getindex(lvw, _kws_to_mask(lattice(lvw), kw))
Base.getindex(l::Lattice; kw...) = getindex(l, _kws_to_mask(l, kw))
Base.setindex!(lvw::LatticeValueWrapper, rhs; kw...) = setindex!x(lvw, rhs, _kws_to_mask(lattice(lvw), kw))

"""
    project(lv::LatticeValue, axis)

Creates a mapping from site coordinates to values of `lv`.
The coordinate axis to project the sites onto can be set with the `axis` argument -
it can be either an integer from 1 to 3 or an axis descriptor `Symbol`.
"""
function project(lv::LatticeValue, axis)
    type, axis_no = parse_axis_sym(axis)
    l = lattice(lv)
    axis_no > dims(l) && error(
        "axis number $axis_no (descriptor '$axis') exceeds lattice dimensionality $(dims(l))")
    crds = collect_coords(l)
    pr_crds = type === Coord ? pr_crds = crds[axis_no, :] :
        crds' * normalize(bravais(l).translation_vectors[:, axis_no])
    perm = sortperm(pr_crds)
    pr_crds[perm], lv.values[perm]
end

ketstate(lv::LatticeValue) = Ket(LatticeBasis(lattice(lv)), lv.values)
brastate(lv::LatticeValue) = Bra(LatticeBasis(lattice(lv)), lv.values)
QuantumOpticsBase.tensor(ket::Ket, lv::LatticeValue) = ket ⊗ ketstate(lv)
QuantumOpticsBase.tensor(bra::Bra, lv::LatticeValue) = bra ⊗ brastate(lv)
QuantumOpticsBase.tensor(lv::LatticeValue, state::QuantumOpticsBase.StateVector) = state ⊗ lv
