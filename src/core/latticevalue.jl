using LinearAlgebra, Logging, QuantumOpticsBase

struct LatticeValueWrapper{LT<:AbstractLattice, VT<:AbstractVector}
    latt::LT
    values::VT
    function LatticeValueWrapper(latt::LT, values::VT) where {LT,VT}
        length(latt) != length(values) &&
            throw(DimensionMismatch("vector has length $(length(values)), lattice has length $(length(latt))"))
        new{LT,VT}(latt, values)
    end
end

lattice(lvw::LatticeValueWrapper) = lvw.latt

Base.copy(lvw::LatticeValueWrapper) = LatticeValueWrapper(lvw.latt, copy(lvw.values))
Base.length(lvw::LatticeValueWrapper) = length(lvw.values)
Base.size(lvw::LatticeValueWrapper) = size(lvw.values)
function Base.getindex(lvw::LatticeValueWrapper, site::AbstractSite)
    i = site_index(lattice(lvw), site)
    i === nothing && throw(BoundsError(lvw, site))
    return lvw.values[i]
end
function Base.setindex!(lvw::LatticeValueWrapper, rhs, site::AbstractSite)
    i = site_index(lattice(lvw), site)
    i === nothing && throw(BoundsError(lvw, site))
    lvw.values[i] = rhs
end
Base.eltype(lvw::LatticeValueWrapper) = eltype(lvw.values)
Base.eachindex(lvw::LatticeValueWrapper) = lattice(lvw)
Base.iterate(lvw::LatticeValueWrapper, s...) = iterate(lvw.values, s...)
Base.pairs(lvw::LatticeValueWrapper) = Iterators.map(=>, lvw.latt, lvw.values)

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
LatticeValue(l::AbstractLattice, v::AbstractVector) = LatticeValueWrapper(l, convert(Vector, v))
LatticeValue(lf, l::AbstractLattice) = LatticeValueWrapper(l, [lf(site) for site in l])

"""
    coord_values(l::Lattice)

Generates a tuple of `LatticeValue`s representing spatial coordinates.
"""
coord_values(l::AbstractLattice) =
    [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]
macro p_str(e::String)
    e === "index" && return :(UnitcellIndex())
    e === "x" && return :(Coord(1))
    e === "y" && return :(Coord(2))
    e === "z" && return :(Coord(3))
    length(e) < 2 && error("Cannot parse site parameter `$e`")
    typ = e[1]
    ind = tryparse(Int, e[2:end])
    ind === nothing && error("Cannot parse site parameter `$e`")
    if typ == 'x'
        return :(Coord($ind))
    elseif typ == 'j'
        return :(UnitcellAxis($ind))
    end
    error("Cannot parse site parameter `$e`")
end
param_value(l::AbstractLattice, a::SiteParameter) = LatticeValue(l, [get_param(site, a) for site in l])
param_value(l::AbstractLattice, a::Symbol) = param_value(l, SiteParameter(a))

Base.rand(l::AbstractLattice) = LatticeValue(l, rand(length(l)))
Base.rand(T::Type, l::AbstractLattice) = LatticeValue(l, rand(T, length(l)))
Base.randn(l::AbstractLattice) = LatticeValue(l, randn(length(l)))
Base.randn(T::Type, l::AbstractLattice) = LatticeValue(l, randn(T, length(l)))
Base.fill(value, l::AbstractLattice) = LatticeValue(l, fill(value, length(l)))
Base.fill!(lv::LatticeValue, value) = (fill!(lv.values, value); lv)
Base.zero(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), zero(lvw.values))
Base.zeros(T::Type, l::AbstractLattice) = fill(zero(T), l)
Base.zeros(l::AbstractLattice) = zeros(Float64, l)
Base.one(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), one(lvw.values))
Base.ones(T::Type, l::AbstractLattice) = fill(one(T), l)
Base.ones(l::AbstractLattice) = ones(Float64, l)

Base.:(==)(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper) = (lvw1.latt == lvw2.latt) && (lvw1.values == lvw2.values)

struct LatticeStyle <: Broadcast.BroadcastStyle end
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{LatticeStyle}) = (copyto!(lvw.values, src); return lvw)
Base.copyto!(lvw::LatticeValueWrapper, src::Broadcast.Broadcasted{Broadcast.DefaultArrayStyle{0}}) = (copyto!(lvw.values, src); return lvw)
Base.setindex!(lvw::LatticeValueWrapper, rhs, i::Int) = setindex!(lvw.values, rhs, i)
Base.broadcastable(lvw::LatticeValueWrapper) = lvw
Base.broadcastable(l::AbstractLattice) = l
Base.getindex(lvw::LatticeValueWrapper, i::CartesianIndex{1}) = lvw.values[only(Tuple(i))]
Base.getindex(l::AbstractLattice, i::CartesianIndex{1}) = l[only(Tuple(i))]
Base.BroadcastStyle(::Type{<:LatticeValueWrapper}) = LatticeStyle()
Base.BroadcastStyle(::Type{<:AbstractLattice}) = LatticeStyle()
Base.BroadcastStyle(bs::Broadcast.BroadcastStyle, ::LatticeStyle) =
    error("cannot broadcast LatticeValue along style $bs")
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LatticeStyle) = LatticeStyle()

function Base.similar(bc::Broadcast.Broadcasted{LatticeStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(bc::Broadcast.Broadcasted) = _extract_lattice(bc.args)
_extract_lattice(ext::Broadcast.Extruded) = _extract_lattice(ext.x)
_extract_lattice(lv::LatticeValueWrapper) = lv.latt
_extract_lattice(x) = x
_extract_lattice(::Tuple{}) = nothing
_extract_lattice(args::Tuple) =
    _extract_lattice(_extract_lattice(args[begin]), Base.tail(args))
_extract_lattice(::Any, rem_args::Tuple) = _extract_lattice(rem_args)
_extract_lattice(l::AbstractLattice, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)

_extract_lattice_s(l::AbstractLattice, args::Tuple) =
    _extract_lattice_s(l, _extract_lattice(args[begin]), Base.tail(args))
_extract_lattice_s(l::AbstractLattice, ::Tuple{}) = l
_extract_lattice_s(l::AbstractLattice, ::Any, rem_args::Tuple) =
    _extract_lattice_s(l, rem_args)
function _extract_lattice_s(l::AbstractLattice, l2::AbstractLattice, rem_args::Tuple)
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

site_inds(l::AbstractLattice{SiteT}, lv_mask::LatticeValue{Bool, <:AbstractLattice{SiteT}}) where SiteT<:AbstractSite =
    Int[site_index(l, lattice(lv_mask)[i])
        for i in eachindex(lv_mask.values) if lv_mask.values[i]]
Base.@propagate_inbounds function Base.getindex(l::AbstractLattice, lv_mask::LatticeValue{Bool})
    @boundscheck check_issublattice(lattice(lv_mask), l)
    inds = site_inds(l, lv_mask)
    l[inds]
end

Base.@propagate_inbounds function Base.getindex(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck check_samelattice(lv, lv_mask)
    inds = site_inds(lattice(lv), lv_mask)
    LatticeValueWrapper(lattice(lv)[inds], lv.values[inds])
end

Base.@propagate_inbounds function Base.Broadcast.dotview(lv::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck check_samelattice(lv, lv_mask)
    inds = site_inds(lattice(lv), lv_mask)
    LatticeValueWrapper(lattice(lv)[inds], @view lv.values[inds])
end
function Base.Broadcast.dotview(lv::LatticeValueWrapper, pairs::Pair{<:SiteParameter}...)
    inds = pairs_to_inds(lattice(lv), pairs...)
    return LatticeValueWrapper(lattice(lv)[inds], @view lv.values[inds])
end
function Base.Broadcast.dotview(lv::LatticeValueWrapper; kw...)
    return Base.Broadcast.dotview(lv, (SiteParameter(crd) => val for (crd, val) in kw)...)
end

Base.@propagate_inbounds function Base.setindex!(lv::LatticeValueWrapper, lv_rhs::LatticeValueWrapper, lv_mask::LatticeValue{Bool})
    @boundscheck begin
        check_issublattice(lattice(lv_mask), lattice(lv))
        check_issublattice(lattice(lv_rhs), lattice(lv))
    end
    inds_l = site_inds(lattice(lv), lv_mask)
    inds_r = site_inds(lattice(lv_rhs), lv_mask)
    lv.values[inds_l] = lv_rhs.values[inds_r]
    return lv_rhs
end

function Base.getindex(lvw::LatticeValueWrapper, pairs::Pair{<:SiteParameter}...)
    inds = pairs_to_inds(lattice(lvw), pairs...)
    return length(inds) == 1 ? lvw.values[only(inds)] :
        LatticeValueWrapper(lattice(lvw)[inds], lvw.values[inds])
end
function Base.getindex(lvw::LatticeValueWrapper; kw...)
    lvw[(SiteParameter(crd) => val for (crd, val) in kw)...]
end
function Base.setindex!(lvw::LatticeValueWrapper, rhs, pairs::Pair{<:SiteParameter}...)
    inds = pairs_to_inds(lattice(lvw), pairs...)
    setindex!(lvw.values, rhs, inds)
end
function Base.setindex!(lvw::LatticeValueWrapper, rhs; kw...)
    lvw[(SiteParameter(crd) => val for (crd, val) in kw)...] = rhs
end

"""
    project(lv::LatticeValue, axis)

Creates a mapping from site coordinates to values of `lv`.
The coordinate axis to project the sites onto can be set with the `axis` argument.
"""
function project(lv::LatticeValue, param::SiteParameter)
    pr_crds = [get_param(site, param) for site in lattice(lv)]
    perm = sortperm(pr_crds)
    pr_crds[perm], lv.values[perm]
end
project(any, param::Symbol) = project(any, SiteParameter(param))

ketstate(lv::LatticeValue) = Ket(LatticeBasis(sites(lv)), lv.values)
brastate(lv::LatticeValue) = Bra(LatticeBasis(sites(lv)), lv.values)
QuantumOpticsBase.tensor(ket::Ket, lv::LatticeValue) = ket ⊗ ketstate(lv)
QuantumOpticsBase.tensor(bra::Bra, lv::LatticeValue) = bra ⊗ brastate(lv)
QuantumOpticsBase.tensor(lv::LatticeValue, state::QuantumOpticsBase.StateVector) = state ⊗ lv
