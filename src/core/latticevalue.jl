using LinearAlgebra, Logging, QuantumOpticsBase

struct LatticeValueWrapper{LT<:AbstractLattice, VT<:AbstractVector}
    lat::LT
    values::VT
    function LatticeValueWrapper(lat::LT, values::VT) where {LT,VT}
        length(lat) != length(values) &&
            throw(DimensionMismatch("vector has length $(length(values)), lattice has length $(length(lat))"))
        new{LT,VT}(lat, values)
    end
end

lattice(lvw::LatticeValueWrapper) = lvw.lat

Base.copy(lvw::LatticeValueWrapper) = LatticeValueWrapper(lvw.lat, copy(lvw.values))
Base.length(lvw::LatticeValueWrapper) = length(lvw.values)
Base.size(lvw::LatticeValueWrapper) = size(lvw.values)
Base.eltype(lvw::LatticeValueWrapper) = eltype(lvw.values)
Base.eachindex(lvw::LatticeValueWrapper) = lattice(lvw)
Base.iterate(lvw::LatticeValueWrapper, s...) = iterate(lvw.values, s...)
Base.pairs(lvw::LatticeValueWrapper) = Iterators.map(=>, lvw.lat, lvw.values)
Base.keys(lvw::LatticeValueWrapper) = lvw.lat
Base.values(lvw::LatticeValueWrapper) = lvw.values

function Base.:(+)(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper)
    check_samelattice(lvw1.lat, lvw2.lat)
    LatticeValueWrapper(lvw1.lat, lvw1.values + lvw2.values)
end
Base.:(-)(lvw1::LatticeValueWrapper) = LatticeValueWrapper(lvw1.lat, -lvw1.values)
function Base.:(-)(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper)
    check_samelattice(lvw1.lat, lvw2.lat)
    LatticeValueWrapper(lvw1.lat, lvw1.values - lvw2.values)
end
Base.:(*)(lvw::LatticeValueWrapper, x::Number) = LatticeValueWrapper(lvw.lat, lvw.values * x)
Base.:(*)(x::Number, lvw::LatticeValueWrapper) = lvw * x
Base.:(/)(lvw::LatticeValueWrapper, x::Number) = LatticeValueWrapper(lvw.lat, lvw.values / x)
Base.isapprox(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper; kw...) =
    isapprox(lvw1.values, lvw2.values; kw...)

"""
    LatticeValue{T, LT}

Represents a value of type `T` on a `LT` lattice.

## Fields
- lattice: the `AbstractLattice` object the value is defined on
- values: the values on different sites
"""
const LatticeValue{T, LT} = LatticeValueWrapper{LT, Vector{T}}

"""
    LatticeValue(lat, values)

Constructs a LatticeValue object.

## Arguments
- `lat`: the lattice the value is defined on.
- `values`: an `AbstractVector` of values on the lattice.

## Example
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(2, 2);

julia> LatticeValue(l, [1, 2, 3, 4])    # Custom values
LatticeValue{Int64} on a 4-site SquareLattice in 2D space
Values stored in a Vector{Int64}:
[1, 2, 3, 4]

julia> x = LatticeValue(l, :x)          # x-coordinate
LatticeValue{Float64} on a 4-site SquareLattice in 2D space
Values stored in a Vector{Float64}:
[1.0, 1.0, 2.0, 2.0]

julia> l2 = SquareLattice(3, 3);

julia> y = LatticeValue(l2, :y)          # y-coordinate
LatticeValue{Float64} on a 9-site SquareLattice in 2D space
Values stored in a Vector{Float64}:
[1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0]

julia> x + y
ERROR: Matching lattices expected.
Got following:
  #1: 4-site SquareLattice in 2D space
  #2: 9-site SquareLattice in 2D space
[...]
```
"""
LatticeValue(l::AbstractLattice, v::AbstractVector) = LatticeValueWrapper(l, convert(Vector, v))
LatticeValue(lf, l::AbstractLattice) = LatticeValueWrapper(l, [lf(site) for site in l])

"""
    coordvalues(lat)

Generates a tuple of `LatticeValue`s representing spatial coordinates on lattice `lat`.
"""
coordvalues(l::AbstractLattice) =
    [LatticeValue(l, vec) for vec in eachrow(collect_coords(l))]

"""
    LatticeValue(lat, prop)

Generates a `LatticeValue` representing the site property `prop` on lattice `lat`.

## Arguments
- `lat`: the lattice the value is defined on.
- `prop`: the `SiteProperty` to be represented. Can be a `SiteProperty` or a `Symbol` defining it.
"""
LatticeValue(l::AbstractLattice, a::SiteProperty) = LatticeValue(l, [getsiteproperty(site, a) for site in l])
LatticeValue(l::AbstractLattice, sym::Symbol) = LatticeValue(l, SitePropertyAlias{sym}())

"""
    coordvalue(lat, coord)

Generates a `LatticeValue` representing the spatial coordinate `coord` on lattice `lat`.

## Arguments
- `lat`: the lattice the value is defined on.
- `coord`: the coordinate to be represented. Can be an integer or a symbol.
"""
coordvalue(l::AbstractLattice, i::Int) = LatticeValue(l, Coord(i))
coordvalue(l::AbstractLattice, sym::Symbol) = LatticeValue(l, SitePropertyAlias{sym}())

Base.rand(l::AbstractLattice) = LatticeValue(l, rand(length(l)))
Base.rand(T::Type, l::AbstractLattice) = LatticeValue(l, rand(T, length(l)))
Base.randn(l::AbstractLattice) = LatticeValue(l, randn(length(l)))
Base.randn(T::Type, l::AbstractLattice) = LatticeValue(l, randn(T, length(l)))
Base.fill(value, l::AbstractLattice) = LatticeValue(l, fill(value, length(l)))
Base.fill!(lv::LatticeValue, value) = (fill!(lv.values, value); lv)
Base.zero(lvw::LatticeValueWrapper) = LatticeValueWrapper(lattice(lvw), zero(lvw.values))
Base.zeros(T::Type, l::AbstractLattice) = fill(zero(T), l)
Base.zeros(l::AbstractLattice) = zeros(Float64, l)
Base.ones(T::Type, l::AbstractLattice) = fill(one(T), l)
Base.ones(l::AbstractLattice) = ones(Float64, l)

Base.:(==)(lvw1::LatticeValueWrapper, lvw2::LatticeValueWrapper) = (lvw1.lat == lvw2.lat) && (lvw1.values == lvw2.values)

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
    throw(ArgumentError("cannot broadcast LatticeValue along $bs"))
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::LatticeStyle) = LatticeStyle()

function Base.similar(bc::Broadcast.Broadcasted{LatticeStyle}, ::Type{Eltype}) where {Eltype}
    l = _extract_lattice(bc)
    LatticeValue(l, similar(Vector{Eltype}, axes(bc)))
end
_extract_lattice(bc::Broadcast.Broadcasted) = _extract_lattice(bc.args)
_extract_lattice(ext::Broadcast.Extruded) = _extract_lattice(ext.x)
_extract_lattice(lv::LatticeValueWrapper) = lv.lat
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

function Base.show(io::IO, ::MIME"text/plain", lv::LatticeValueWrapper)
    print(io, "LatticeValue{$(eltype(lv))} on a ")
    summary(io, lattice(lv))
    if !requires_compact(io)
        print(io, "\nValues stored in a $(typeof(lv.values)):\n")
        Base.show_vector(IOContext(io, :compact => true), lv.values)
    end
end

function to_inds(l::AbstractLattice, sites::AbstractVector{<:AbstractSite})
    indices = Int[]
    for site in sites
        index = site_index(l, site)
        index === nothing && return nothing
        push!(indices, index)
    end
    return indices
end
function to_inds(l::AbstractLattice, lv_mask::LatticeValue{Bool})
    indices = Int[]
    l2 = lattice(lv_mask)
    check_issublattice(l, l2)
    for (i, site) in enumerate(l2)
        if lv_mask.values[i]
            index = site_index(l, site)
            index === nothing && continue
            push!(indices, index)
        end
    end
    return indices
end
function to_inds(l::AbstractLattice, l2::AbstractLattice)
    check_issublattice(l2, l)
    Int[site_index(l, site) for site in l2]
end
function to_inds(l::AbstractLattice, pairs::Pair...; kw...)
    inds = pairs_to_indices(l, to_param_pairs(pairs...; kw...))
    checkbounds(l, inds)
    return inds
end
function to_inds(l::AbstractLattice, ::typeof(!), pairs::Pair...; kw...)
    ind = pairs_to_index(l, to_param_pairs(pairs...; kw...))
    return ind
end
to_inds(l::AbstractLattice, site::AbstractSite) = site_index(l, site)
to_inds(::AbstractLattice, rs::ResolvedSite) = rs.index

Base.@propagate_inbounds function Base.getindex(l::AbstractLattice, args...; kw...)
    inds = to_inds(l, args...; kw...)
    if inds isa Nothing || isempty(inds)
        throw(BoundsError(l, (args..., NamedTuple(kw))))
    elseif length(inds) == 1
        return l[only(inds)]
    else
        return l[inds]
    end
end
Base.@propagate_inbounds function Base.getindex(lv::LatticeValueWrapper, args...; kw...)
    inds = to_inds(lattice(lv), args...; kw...)
    if inds isa Nothing || isempty(inds)
        throw(BoundsError(lv, (args..., NamedTuple(kw))))
    elseif length(inds) == 1
        return lv.values[only(inds)]
    else
        return LatticeValueWrapper(lattice(lv)[inds], lv.values[inds])
    end
end

Base.@propagate_inbounds function Base.view(lv::LatticeValueWrapper, args...; kw...)
    inds = to_inds(lattice(lv), args...; kw...)
    if inds === nothing
        throw(BoundsError(lv, (args..., NamedTuple(kw))))
    end
    LatticeValueWrapper(lattice(lv)[inds], @view lv.values[inds])
end

Base.@propagate_inbounds Base.Broadcast.dotview(lv::LatticeValueWrapper, args...; kw...) =
    view(lv, args...; kw...)

Base.@propagate_inbounds function Base.setindex!(lv::LatticeValueWrapper, lv_rhs::LatticeValueWrapper, args...; kw...)
    inds_l = to_inds(lattice(lv), args...; kw...)
    inds_l === nothing && throw(BoundsError(lv, (args..., NamedTuple(kw))))
    inds_r = Int[]
    for i in inds_l
        index = site_index(lattice(lv_rhs), lattice(lv)[i])
        index === nothing && throw(ArgumentError("Cannot assign: site not found in lattice"))
        push!(inds_r, index)
    end
    lv.values[inds_l] = @view lv_rhs.values[inds_r]
    return lv_rhs
end
function Base.setindex!(lvw::LatticeValueWrapper, rhs::Any, args...; kw...)
    inds = to_inds(lattice(lvw), args...; kw...)
    if inds === nothing
        throw(BoundsError(lvw, (args..., NamedTuple(kw))))
    elseif inds isa Int
        lvw.values[inds] = rhs
    else
        # throw an error, because we don't want to broadcast
        lvw.values[inds] = 0
    end
    return rhs
end

"""
    ketstate(lv)

Converts a `LatticeValue` to a `Ket` wavefunction vector.
"""
ketstate(lv::LatticeValue) = Ket(LatticeBasis(lattice(lv)), lv.values)

"""
    brastate(lv)

Converts a `LatticeValue` to a `Bra` wavefunction vector.
"""
brastate(lv::LatticeValue) = Bra(LatticeBasis(lattice(lv)), lv.values)
QuantumOpticsBase.tensor(ket::Ket, lv::LatticeValue) = ket ⊗ ketstate(lv)
QuantumOpticsBase.tensor(bra::Bra, lv::LatticeValue) = bra ⊗ brastate(lv)
QuantumOpticsBase.tensor(lv::LatticeValue, state::QuantumOpticsBase.StateVector) = state ⊗ lv
