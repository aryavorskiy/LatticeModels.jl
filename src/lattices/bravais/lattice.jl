using LinearAlgebra, Logging, StaticArrays

struct BravaisLattice{N, UnitcellT, BS<:BoundaryConditions} <:
        AbstractLattice{BravaisSite{N, UnitcellT}}
    unitcell::UnitcellT
    pointers::Vector{BravaisPointer{N}}
    boundaries::BS
    b_depth::Int
    function BravaisLattice(unitcell::UnitcellT, pointers::Vector{BravaisPointer{N}}, boundaries::BoundariesT;
            b_depth=1) where {N,UnitcellT,BoundariesT}
        dims(unitcell) != N && throw(ArgumentError("Dimension mismatch"))
        new{N,UnitcellT,BoundariesT}(unitcell, sort!(pointers), boundaries, b_depth)
    end
end
BravaisLattice(unitcell, pointers) = BravaisLattice(unitcell, pointers, BoundaryConditions())
add_boundaries(l::BravaisLattice, bs) = BravaisLattice(l.unitcell, l.pointers, to_boundaries(bs))
add_boundaries(l::BravaisLattice, ::Nothing) = l
add_boundaries(l::AbstractLattice, ::Nothing) = l
sites(l::BravaisLattice) = Sites(add_boundaries(l, BoundaryConditions()))
sublatvector(l::BravaisLattice, i::Int) = sublatvector(l.unitcell, i)
trvectors(l::BravaisLattice) = trvectors(l.unitcell)
b_depth(l::BravaisLattice) = l.b_depth
cartesian_indices(l::BravaisLattice{N, B, <:BoundaryConditions{<:NTuple{M}}} where {N, B}) where M =
    cartesian_indices(l.b_depth, Val(M))

Base.:(==)(l1::BravaisLattice, l2::BravaisLattice) =
    (l1.pointers == l2.pointers) && (l1.unitcell == l2.unitcell)

# Set functions
Base.emptymutable(l::BravaisLattice{N, B}, ::Type{BravaisSite{N, B}}) where {N, B} =
    BravaisLattice(l.unitcell, BravaisPointer{N}[])
Base.copymutable(l::BravaisLattice) = BravaisLattice(l.unitcell, copy(l.pointers))
Base.length(l::BravaisLattice) = length(l.pointers)
basis_length(l::BravaisLattice) = length(l.unitcell)
function check_sameunitcell(l::BravaisLattice, site::BravaisSite)
    l.unitcell != site.unitcell && throw(ArgumentError("Unit cell mismatch: $(l.unitcell) ≠ $(site.unitcell)"))
end

get_site(l::BravaisLattice, lp::BravaisPointer) = BravaisSite(lp, l.unitcell)
get_site(::BravaisLattice, site::BravaisSite) = site
get_site(::BravaisLattice, ::NoSite) = NoSite()

function resolve_site(l::BravaisLattice, site::AbstractSite)
    factor, site = shift_site(l.boundaries, l, site)
    i = site_index(l, site)
    i === nothing && return nothing
    return ResolvedSite(site, i, factor)
end

Base.in(lp::BravaisPointer, l::BravaisLattice) = insorted(lp, l.pointers)
function Base.in(site::BravaisSite{N, B}, l::BravaisLattice{N, B}) where {N, B}
    @boundscheck check_sameunitcell(l, site)
    return in(site.lp, l)
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
    @assert 1 ≤ lp.basis_index ≤ basis_length(l) "invalid basis index $(lp.basis_index)"
    i = searchsortedfirst(l.pointers, lp)
    i ≤ length(l) && l.pointers[i] == lp && return l
    insert!(l.pointers, i, lp)
end
Base.push!(l::BravaisLattice{N, B}, site::BravaisSite{N, B}) where {N, B} =
    push!(l, site.lp)

Base.@propagate_inbounds function site_index(l::BravaisLattice, lp::BravaisPointer)
    i = searchsortedfirst(l.pointers, lp)
    i > length(l) && return nothing
    return l.pointers[i] == lp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::BravaisLattice, site::BravaisSite)
    @boundscheck check_sameunitcell(l, site)
    site_index(l, site.lp)
end

function Base.summary(io::IO, l::BravaisLattice{N, <:UnitCell{Sym}}) where {N,Sym}
    print(io, length(l), "-site ", N, "-dim ", Sym)
    if basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", l::BravaisLattice)
    summary(io, l)
    (isempty(l.boundaries.bcs) || get(io, :compact, false)) && return
    print(io, " with boundary conditions:")
    for bc in l.boundaries.bcs
        println()
        show(io, mime, bc)
    end
end

const RangeT = Union{Integer, OrdinalRange{<:Integer, <:Integer}}
function _sort(i::Integer)
    i > 0 && return i
    throw(ArgumentError("Invalid range: $i; must be positive"))
end
function _sort(r::OrdinalRange)
    if step(r) > 0
        return r
    else
        return reverse(r)
    end
end

"""
    span_unitcells(unitcell, dims...[; kw...])

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

julia> uc = UnitCell([[1, 0]; [0, 1]]);

julia> span_unitcells(uc, 3, 3) == SquareLattice(3, 3)
true
```
"""
function span_unitcells(unitcell::UnitCell{Sym,N,NB}, sz::Vararg{RangeT, N};
        boundaries=BoundaryConditions(), offset = :origin) where {Sym,N,NB}
    ptrs = BravaisPointer{N}[]
    for unit_cell in CartesianIndices(reverse(_sort.(sz)))
        svec = reverse(SVector{N}(Tuple(unit_cell)))
        for i in 1:NB
            push!(ptrs, BravaisPointer(svec, i))
        end
    end
    return BravaisLattice(offset_unitcell(unitcell, offset), ptrs, to_boundaries(boundaries))
end
span_unitcells(uc::UnitCell, sz::Tuple{Vararg{RangeT}}; kw...) =
    throw(ArgumentError("Dimension mismatch: $(dims(uc))-dim unit cell, $(length(sz)) lattice dimensions"))
