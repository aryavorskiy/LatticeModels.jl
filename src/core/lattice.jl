using LinearAlgebra, Logging, StaticArrays

"""
    Lattice{N,B}
A finite subset of a `B` bravais lattice.

---
    Lattice(sym, sz, bvs[, mask])
Constructs a finite `Lattice{sym, N, NB}` as a subset of the `bvs` Bravais lattice.
`sz` is a `NTuple{N, Int}` which represents how many times the unit cell of `bvs` was translated by each axis - these sites form a *macrocell*.
`mask`, if defined, is a `Vector{Bool}` storing information about which of the sites from the macrocell
are actually included in the lattice, and which are not.

For example, a 3×3 square lattice with its center site excluded is represented as
`Lattice(:square, (3, 3), Bravais([1 0; 0 1]), Bool[1, 1, 1, 1, 0, 1, 1, 1, 1])`

To define a new type of lattice, create an alias for `Lattice{YourSym, YourN, YourNB}`.
Refer to the docs for detailed explanation.
"""
struct BravaisLattice{N, B<:UnitCell{Sym, N} where Sym, BS} <: AbstractLattice{BravaisSite{N, B}}
    bravais::B
    pointers::Vector{BravaisPointer{N}}
    boundaries::BS
    function BravaisLattice(bravais::B, pointers::Vector{BravaisPointer{N}}, boundaries::BS) where {N,B,BS}
        new{N,B,BS}(bravais, sort!(pointers), boundaries)
    end
end
BravaisLattice(bravais, pointers) = BravaisLattice(bravais, pointers, ())

Base.:(==)(l1::BravaisLattice, l2::BravaisLattice) = (l1.pointers == l2.pointers) && (l1.bravais == l2.bravais)

Base.emptymutable(l::BravaisLattice{B, N}, ::Type{BravaisSite{N}}=eltype(l)) where {B, N} =
    BravaisLattice(l.bravais, [])
Base.copymutable(l::BravaisLattice) = BravaisLattice(l.bravais, copy(l.pointers))
Base.length(l::BravaisLattice) = length(l.pointers)
lattice_type(::BravaisLattice{<:UnitCell{Sym}}) where {Sym} = Sym
basis_length(l::BravaisLattice) = length(l.bravais)
check_bravais(l::BravaisLattice, site::BravaisSite) =
    @assert l.bravais == site.bravais

default_bonds(::BravaisLattice, ::Val) = ()
default_bonds(l::BravaisLattice) = default_bonds(l, Val(1))
default_bonds(l::BravaisLattice, i::Int) = default_bonds(l, Val(i))

get_site(l::BravaisLattice, lp::BravaisPointer) = BravaisSite(lp, l.bravais)
get_site(::BravaisLattice, site::BravaisSite) = site
get_site(::BravaisLattice, ::Nothing) = nothing

Base.in(l::BravaisLattice, lp::BravaisPointer) = insorted(lp, l.pointers)
function Base.in(l::BravaisLattice{N, B}, site::BravaisSite{N, B}) where {N, B}
    @boundscheck check_samebravais(l, site)
    in(l, site.lp)
end

Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, i::Int)
    @boundscheck checkbounds(l, i)
    return get_site(l, l.pointers[i])
end
Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, is::AbstractVector{Int})
    @boundscheck checkbounds(l, is)
    return BravaisLattice(l.bravais, l.pointers[sort(is)])
end
function Base.delete!(l::BravaisLattice, lp::BravaisPointer)
    i = site_index(l, lp)
    i !== nothing && deleteat!(l.pointers, i)
    l
end
function Base.delete!(l::BravaisLattice{N, B}, site::BravaisSite{N, B}) where {N, B}
    check_samebravais(l, site)
    Base.delete!(l, site.lp)
end
function Base.deleteat!(l::BravaisLattice, inds)
    @boundscheck checkbounds(l, inds)
    deleteat!(l.pointers, inds)
    l
end

function Base.push!(l::BravaisLattice, lp::BravaisPointer)
    @assert 1 ≤ lp.basis_index ≤ basis_length(lp) "invalid basis index $(lp.basis_index)"
    i = searchsortedfirst(l.pointers, lp)
    i < length(l) && l.pointers[i] != lp && insert(l.pointers, i, lp)
end
Base.push!(l::BravaisLattice{N, B}, site::BravaisSite{N, B}) where {N, B} =
    push!(l, site.lp)

"""
    site_index(l::Lattice, site::LatticeSite; macrocell=false)

Returns the integer index for given `site` in `lattice`.
Returns `nothing` if the site is not present in the lattice.
"""
Base.@propagate_inbounds function site_index(l::BravaisLattice, lp::BravaisPointer)
    i = searchsortedfirst(l.pointers, lp)
    i > length(l) && return nothing
    return l.pointers[i] == lp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::BravaisLattice, site::BravaisSite)
    @boundscheck l.bravais == site.bravais
    site_index(l, site.lp)
end

function Base.iterate(l::BravaisLattice, state = (1, length(l)))
    i, len = state
    return i > len ? nothing : (l[i], (i+1, len))
end

function Base.show(io::IO, ::MIME"text/plain", l::BravaisLattice{N, <:UnitCell{Sym}}) where {N,Sym}
    print(io, length(l), "-site ", N, "-dimensional ", Sym, " lattice")
    if basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

function BravaisLattice(bvs::UnitCell{Sym,N,NB}, sz::NTuple{N,Int}) where {Sym,N,NB}
    ptrs = BravaisPointer{N}[]
    for unit_cell in CartesianIndex{N}():CartesianIndex(reverse(sz))
        svec = reverse(SVector{N}(Tuple(unit_cell)))
        for i in 1:NB
            push!(ptrs, BravaisPointer(svec, i))
        end
    end
    BravaisLattice(bvs, ptrs)
end

const InfDimLattice{Sym,N,NB} = BravaisLattice{N, <:UnitCell{Sym,N,NB}}
const AnyDimLattice{Sym,NB} = BravaisLattice{N, <:UnitCell{Sym,N,NB}} where N
function InfDimLattice{Sym,N,NB}(sz::Vararg{Int,N}) where {Sym,N,NB}
    return BravaisLattice(UnitCell{Sym,N,NB}(), sz)
end
function AnyDimLattice{Sym,NB}(sz::Vararg{Int,N}) where {Sym,N,NB}
    return BravaisLattice(UnitCell{Sym,N,NB}(), sz)
end
(::Type{T})(f::Function, sz::Vararg{Int}) where T<:BravaisLattice =
    sublattice(f, T(sz...))

macrocell_size(::BravaisLattice) = error("This function is discontinued")
