using Logging, StaticArrays

"""
    UnitCell{Sym, N, NB}
`N`-dimensional infinite Bravais lattice unit cell with `NB` sites in basis.
`Sym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Bravais` object behavior known at compile-time,
which allows to introduce various optimizations or define specific plot recipes.

---
    UnitCell(translation_vectors[, basis, offset])
Constructs a Bravais lattice unit cell with given translation vectors and locations of basis sites
relative to some unit cell.
The `basis` argument can be omitted, in which case the lattice basis will consist of
one site located in the bottom-left corner of the unit cell.

`translation_vectors` argument must be an `AbstractMatrix{<:Real}` of size `N×N`,
while `basis` must also be an  abstract matrix of size `N×NB`.
"""
struct UnitCell{Sym,N,NB,NN,NNB}
    translation_vectors::SMatrix{N,N,Float64,NN}
    basis::SMatrix{N,NB,Float64,NNB}
    function UnitCell{Sym}(translation_vectors::AbstractMatrix{<:Real},
        basis::AbstractMatrix{<:Real}=zeros((size(translation_vectors, 1), 1)),
        offset::AbstractVector{<:Real}=zeros(size(translation_vectors, 1))) where Sym
        @check_size translation_vectors :square
        N = size(translation_vectors, 1)
        @check_size offset N
        NB = size(basis, 2)
        @check_size basis (N, NB)
        new{Sym,N,NB,N*N,N*NB}(translation_vectors, basis .+ offset)
    end
end
UnitCell(uc::UnitCell{Sym, N}, offset::SVector{N}) where {Sym, N} =
    UnitCell{Sym}(uc.translation_vectors, uc.basis, offset)

dims(@nospecialize _::UnitCell{Sym, N} where Sym) where {N} = N
Base.:(==)(b1::BT, b2::BT) where BT<:UnitCell =
    b1.translation_vectors == b2.translation_vectors && b1.basis == b2.basis
Base.length(::UnitCell{Sym,N,NB} where {Sym,N}) where {NB} = NB

struct BravaisPointer{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end
dims(::BravaisPointer{N}) where N = N

site_coords(uc::UnitCell, lp::BravaisPointer) =
    uc.basis[:, lp.basis_index] + uc.translation_vectors * lp.unit_cell
site_coords(b::UnitCell{Sym,N,1} where {Sym}, lp::BravaisPointer{N}) where {N} =
    vec(b.basis) + b.translation_vectors * lp.unit_cell

Base.:(==)(lp1::BravaisPointer, lp2::BravaisPointer) =
    lp1.basis_index == lp2.basis_index && lp1.unit_cell == lp2.unit_cell
function Base.isless(site1::BravaisPointer, site2::BravaisPointer)
    if site1.unit_cell == site2.unit_cell
        return isless(site1.basis_index, site2.basis_index)
    else
        return isless(site1.unit_cell, site2.unit_cell)
    end
end

"""
    BravaisSite{N, B}
A site of a `BravaisLattice{N, B}` lattice.

Fields:
- `unit_cell`: a set of translations along all axes representing the unit cell the site is located in.
- `basis_index`: the number of site in the lattice basis.

This type is used to iterate over all sites of a `Lattice{LatticeSym, N, NB}`.
The exact location of a `LatticeSite` can be found using the `site.coords` function.
"""
struct BravaisSite{N, B} <: AbstractSite{N}
    lp::BravaisPointer{N}
    bravais::B
    coords::SVector{N,Float64}
    BravaisSite(lp::BravaisPointer{N}, b::B) where {N, B<:UnitCell} =
        new{N, B}(lp, b, site_coords(b, lp))
end
BravaisSite(::Nothing, ::Any) = NoSite()

@inline function Base.getproperty(site::BravaisSite{N}, sym::Symbol) where N
    if sym === :index
        return site.lp.basis_index
    elseif sym === :unit_cell
        return site.lp.unit_cell
    elseif sym === :x && N ≥ 1
        return site.coords[1]
    elseif sym === :y && N ≥ 2
        return site.coords[2]
    elseif sym === :z && N ≥ 3
        return site.coords[3]
    end
    Base.getfield(site, sym)
end

Base.show(io::IO, ::MIME"text/plain", site::BravaisSite{N, <:UnitCell{Sym}}) where {N, Sym} =
    print(io, "Site of a ", N, "-dim ", Sym, " lattice @ x = $(site.coords)")

struct UnitcellAxis <: SiteParameter axis::Int end
struct UnitcellIndex <: SiteParameter end
function get_param(site::BravaisSite, c::UnitcellAxis)
    @assert 1 ≤ c.axis ≤ dims(site)
    return site.unit_cell[c.axis]
end
get_param(site::BravaisSite, ::UnitcellIndex) = site.lp.basis_index

Base.:(==)(site1::BravaisSite, site2::BravaisSite) =
    site1.lp == site2.lp && site1.coords == site2.coords
Base.isless(site1::BravaisSite, site2::BravaisSite) =
    isless(site1.lp, site2.lp)
