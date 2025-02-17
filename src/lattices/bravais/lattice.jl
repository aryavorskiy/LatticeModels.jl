using LinearAlgebra, Logging, StaticArrays

struct BravaisLattice{N,NU,UnitcellT} <: AbstractLattice{BravaisSite{N,NU,UnitcellT}}
    unitcell::UnitcellT
    pointers::Vector{BravaisPointer{NU}}
    function BravaisLattice(unitcell::UnitcellT, pointers::Vector{BravaisPointer{NU}}) where {N,NU,UnitcellT<:UnitCell{N,_NU} where _NU}
        ldims(unitcell) != NU && throw(ArgumentError("Dimension mismatch: $(ldims(unitcell)) ≠ $NU"))
        !issorted(pointers) && throw(ArgumentError("`pointers` must be sorted"))
        for j in 1:length(pointers) - 1
            pointers[j] == pointers[j + 1] &&
                throw(ArgumentError("Duplicate pointer at index $j: $(pointers[j])"))
        end
        new{N,NU,UnitcellT}(unitcell, pointers)
    end
end
unitcell(l::BravaisLattice) = l.unitcell
unitcell(lw::LatticeWithMetadata) = unitcell(lw.lat)
basvector(any, i::Int) = basvector(unitcell(any), i)
unitvector(any, i::Int) = unitvector(unitcell(any), i)
unitvectors(any) = unitvectors(unitcell(any))
baslength(any) = length(unitcell(any))
transform_unitcell(l::BravaisLattice; kw...) = BravaisLattice(transform_unitcell(l.unitcell; kw...), l.pointers)
transform_unitcell(lw::LatticeWithMetadata{<:BravaisLattice}; kw...) =
    LatticeWithMetadata(transform_unitcell(lw.lat; kw...), lw.metadata)
_showbonds_default(::BravaisLattice) = 1

Base.:(==)(l1::BravaisLattice, l2::BravaisLattice) =
    (l1.pointers == l2.pointers) && (l1.unitcell == l2.unitcell)

# Set functions
Base.emptymutable(l::BravaisLattice{N,NU,B}, ::Type{BravaisSite{N,NU,B}}) where {N,NU,B} =
    BravaisLattice(l.unitcell, BravaisPointer{NU}[])
Base.copymutable(l::BravaisLattice) = BravaisLattice(l.unitcell, copy(l.pointers))
Base.length(l::BravaisLattice) = length(l.pointers)
function check_sameunitcell(l::BravaisLattice, site::BravaisSite)
    l.unitcell != site.unitcell && throw(ArgumentError("Unit cell mismatch: $(l.unitcell) ≠ $(site.unitcell)"))
end

Base.in(bp::BravaisPointer, l::BravaisLattice) = insorted(bp, l.pointers)
function Base.in(site::BravaisSite{N,NU,B}, l::BravaisLattice{N,NU,B}) where {N,NU,B}
    @boundscheck check_sameunitcell(l, site)
    return in(bravaispointer(site), l)
end

Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, i::Int)
    @boundscheck checkbounds(l, i)
    return BravaisSite(l.pointers[i], l.unitcell)
end
Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, is::AbstractVector{Int})
    @boundscheck checkbounds(l, is)
    return BravaisLattice(l.unitcell, l.pointers[sort(is)])
end
Base.@propagate_inbounds function Base.deleteat!(l::BravaisLattice, inds)
    @boundscheck checkbounds(l, inds)
    deleteat!(l.pointers, inds)
    return l
end

function Base.push!(l::BravaisLattice, bp::BravaisPointer)
    @assert 1 ≤ bp.basindex ≤ baslength(l) "invalid basis index $(bp.basindex)"
    i = searchsortedfirst(l.pointers, bp)
    i ≤ length(l) && l.pointers[i] == bp && return l
    insert!(l.pointers, i, bp)
    return l
end
Base.push!(l::BravaisLattice{N,NU,B}, site::BravaisSite{N,NU,B}) where {N,NU,B} =
    push!(l, bravaispointer(site))

Base.@propagate_inbounds function site_index(l::BravaisLattice, bp::BravaisPointer, range)
    ilo, ihi = range[begin], range[end]
    i = searchsortedfirst(l.pointers, bp, ilo, ihi, Base.Order.Forward)
    i in range || return nothing
    return l.pointers[i] == bp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::BravaisLattice, site::BravaisSite, range)
    @boundscheck check_sameunitcell(l, site)
    site_index(l, bravaispointer(site), range)
end

latticename(::BravaisLattice{N,NU,<:UnitCell{N,NU,NB}} where N) where {NU,NB} =
    "$(NU)-dim Bravais lattice" * (NB > 1 ? " with $NB-site basis" : "")
function Base.show(io::IO, mime::MIME"text/plain", l::BravaisLattice)
    summary(io, l)
    if !requires_compact(io)
        println(io, "\nUnit cell:")
        io = IOContext(io, :showtitle => false)
        show(io, mime, l.unitcell)
    end
end

const RangeT = Union{Integer, OrdinalRange{<:Integer, <:Integer}}
function _sort(i::Integer)
    i > 0 && return i
    throw(ArgumentError("Invalid range: $i; must be positive"))
end
function _sort(r::OrdinalRange)
    step(r) < 0 && return reverse(r)
    return r
end

function add_bravaispointers!(f, ptrs, unitcell::UnitCell{N,NU,NB},
        axes::Tuple{Vararg{RangeT}}) where {N,NU,NB}
    for latcoords in CartesianIndices(reverse(_sort.(axes)))
        svec = reverse(SVector{NU}(Tuple(latcoords)))
        for i in 1:NB
            bp = BravaisPointer(svec, i)
            site = BravaisSite(bp, unitcell)
            f(site) && push!(ptrs, BravaisPointer(svec, i))
        end
    end
end

function finalize_lattice(lat; boundaries=BoundaryConditions(), default_translations=(),
        rmdup=false, postoffset=:origin, postrotate=nothing)
    lat = setnnbonds(lat, getnnbonds(stripmeta(lat)))
    lat = addtranslations(lat, default_translations)
    lat = setboundaries(lat, parse_boundaries(lat, boundaries), rmdup=rmdup)
    lat = transform_unitcell(lat, offset=postoffset, rotate=postrotate)
    return lat
end

"""
    span_unitcells([f, ]unitcell, dims...[; boundaries, offset])

Construct a Bravais lattice by spanning `unitcell` in `dims` dimensions, filtered by `f`.
## Arguments
- `f`: a function that defines if the site is included in the lattice. Takes a `BravaisSite`, returns a `Bool`.
- `unitcell`: a `UnitCell` object.
- `dims`: a list of `Integer`s or `Range`s specifying the size of the lattice in each dimension.

## Keyword arguments
- `default_translations`: a list of `BravaisTranslation`s to add to the lattice as default boundary condition axes.
- `boundaries`: a `BoundaryConditions` object specifying the boundary conditions of the lattice.
- `rmdup`: a `Bool` specifying whether to remove sites that are equivalent after applying the boundary conditions.
- `offset`: the offset of the lattice from the origin. See `UnitCell` for details.
- `rotate`: a rotation matrix to apply to the lattice. See `UnitCell` for details.

Keep in mind that the offset and rotation are applied to the unit cell before the lattice is
spanned (and `f` is applied). To apply them after the lattice is spanned, use the `postoffset`
and `postrotate` keywords.

## Examples
```jldoctest
julia> using LatticeModels

julia> using LatticeModels

julia> uc = UnitCell([[1, 0] [0, 1]])
1-site Unit cell of a 2-dim Bravais lattice in 2D space:
  Basis site coordinates:
    ⎡ 0.000⎤
    ⎣ 0.000⎦
  Translation vectors:
    ⎡ 1.000⎤  ⎡ 0.000⎤
    ⎣ 0.000⎦  ⎣ 1.000⎦

julia> span_unitcells(uc, 3, 3) == SquareLattice(3, 3)
true
```
"""
function span_unitcells(f, unitcell::UnitCell{N,NU,NB}, axes::Vararg{RangeT,NU};
        unitvectortrs=true, offset=:origin, rotate=nothing, kw...) where {N,NU,NB}
    ptrs = BravaisPointer{NU}[]
    new_unitcell = transform_unitcell(unitcell, rotate=rotate, offset=offset)
    add_bravaispointers!(f, ptrs, new_unitcell, axes)
    b = BravaisLattice(new_unitcell, ptrs)
    unitvectortrs && for (i, ax) in enumerate(axes)
        nh = ax isa Integer ? ax : length(ax)
        b = addtranslations(b, Symbol("axis$i") => BravaisTranslation(one_hot(i, NU) * nh))
    end
    return finalize_lattice(b; kw...)
end
span_unitcells(f, uc::UnitCell, sz::Vararg{RangeT}; kw...) =
    throw(ArgumentError("Dimension mismatch: $(ldims(uc))-dim unit cell, $(length(sz)) lattice dimensions"))
@inline alwaystrue(x) = true
span_unitcells(uc::UnitCell, axes...; kw...) = span_unitcells(alwaystrue, uc, axes...; kw...)
