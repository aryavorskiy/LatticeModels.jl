using LinearAlgebra, Logging, StaticArrays

struct BravaisLattice{N,NU,UnitcellT} <: AbstractLattice{BravaisSite{N,NU,UnitcellT}}
    unitcell::UnitcellT
    pointers::Vector{BravaisPointer{NU}}
    function BravaisLattice(unitcell::UnitcellT, pointers::Vector{BravaisPointer{NU}}) where {N,NU,UnitcellT<:UnitCell{N,_NU} where _NU}
        ldims(unitcell) != NU && throw(ArgumentError("Dimension mismatch: $(ldims(unitcell)) ≠ $NU"))
        new{N,NU,UnitcellT}(unitcell, sort!(pointers))
    end
end
unitcell(l::BravaisLattice) = l.unitcell
unitcell(lw::LatticeWithParams) = unitcell(lw.lat)
unitvectors(any) = unitvectors(unitcell(any))
basvector(any, i::Int) = basvector(unitcell(any), i)
baslength(any) = length(unitcell(any))
lattransform(ltr::LatticeTransform, l::BravaisLattice) =
    BravaisLattice(l.unitcell |> ltr, unique(lp |> ltr for lp in l.pointers))

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

get_site(l::BravaisLattice, lp::BravaisPointer) = BravaisSite(lp, l.unitcell)
get_site(::BravaisLattice, site::BravaisSite) = site
get_site(::BravaisLattice, ::NoSite) = NoSite()

Base.in(lp::BravaisPointer, l::BravaisLattice) = insorted(lp, l.pointers)
function Base.in(site::BravaisSite{N,NU,B}, l::BravaisLattice{N,NU,B}) where {N,NU,B}
    @boundscheck check_sameunitcell(l, site)
    return in(bravaispointer(site), l)
end

Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, i::Int)
    @boundscheck checkbounds(l, i)
    return get_site(l, l.pointers[i])
end
Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, is::AbstractVector{Int})
    @boundscheck checkbounds(l, is)
    return BravaisLattice(l.unitcell, l.pointers[sort(is)])
end
Base.@propagate_inbounds function Base.delete!(l::BravaisLattice, lp::BravaisPointer)
    i = site_index(l, lp)
    i !== nothing && deleteat!(l.pointers, i)
    return l
end
Base.@propagate_inbounds function Base.deleteat!(l::BravaisLattice, inds)
    @boundscheck checkbounds(l, inds)
    deleteat!(l.pointers, inds)
    return l
end

function Base.push!(l::BravaisLattice, lp::BravaisPointer)
    @assert 1 ≤ lp.basindex ≤ baslength(l) "invalid basis index $(lp.basindex)"
    i = searchsortedfirst(l.pointers, lp)
    i ≤ length(l) && l.pointers[i] == lp && return l
    insert!(l.pointers, i, lp)
    return l
end
Base.push!(l::BravaisLattice{N,NU,B}, site::BravaisSite{N,NU,B}) where {N,NU,B} =
    push!(l, bravaispointer(site))

Base.@propagate_inbounds function site_index(l::BravaisLattice, lp::BravaisPointer, range)
    i = searchsortedfirst(@view(l.pointers[range]), lp) + first(range) - 1
    i in range || return nothing
    return l.pointers[i] == lp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::BravaisLattice, site::BravaisSite, range)
    @boundscheck check_sameunitcell(l, site)
    site_index(l, bravaispointer(site), range)
end

function Base.summary(io::IO, l::BravaisLattice{N,NU,<:UnitCell{N,NU,NB}}) where {N,NU,NB}
    print(io, length(l), "-site $NU-dim Bravais lattice in $(N)D space")
    if baslength(l) > 1
        print(io, " ($NB-site basis)")
    end
end

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
    step(r) > 0 && return reverse(r)
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

function finalize_lattice(lat; boundaries=BoundaryConditions(), default_translations=(), rmdup=false)
    lat = setnnbonds(lat, getnnbonds(stripparams(lat)))
    lat = addtranslations(lat, default_translations)
    lat = setboundaries(lat, parse_boundaries(lat, boundaries), rmdup=rmdup)
    return lat
end

"""
    span_unitcells(unitcell, dims...[; boundaries, offset])

Construct a Bravais lattice by spanning `unitcell` in `dims` dimensions.

## Arguments
- `unitcell`: a `UnitCell` object.
- `dims`: a list of `Integer`s or `Range`s specifying the size of the lattice in each dimension.

## Keyword arguments
- `boundaries`: a `BoundaryConditions` object specifying the boundary conditions of the lattice.
- `offset`: the offset of the lattice from the origin. See `UnitCell` for details.

## Examples
```jldoctest
julia> using LatticeModels

julia> uc = UnitCell([[1, 0] [0, 1]]);

julia> span_unitcells(uc, 3, 3) == SquareLattice(3, 3)
true
```
"""
function span_unitcells(f, unitcell::UnitCell{N,NU,NB}, axes::Vararg{RangeT,NU};
        unitvectortrs=true, offset = :origin, kw...) where {N,NU,NB}
    ptrs = BravaisPointer{NU}[]
    new_unitcell = offset_unitcell(unitcell, offset)
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
