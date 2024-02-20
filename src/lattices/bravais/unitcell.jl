using Logging, StaticArrays, Printf

"""
    UnitCell{Sym,N,NB}
`N`-dimensional infinite Bravais lattice unit cell with `NB` sites in basis.
`Sym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Bravais` object behavior known at compile-time,
which allows to introduce various optimizations or define specific plot recipes.
"""
struct UnitCell{Sym,N,NB,NN,NNB}
    translations::SMatrix{N,N,Float64,NN}
    basissites::SMatrix{N,NB,Float64,NNB}
    function UnitCell{Sym}(translations::AbstractMatrix,
            basissites::AbstractMatrix=zeros(size(translations, 1), 1)) where {Sym}
        @check_size translations :square
        N = size(translations, 1)
        NB = size(basissites, 2)
        @check_size basissites (N, NB)
        new{Sym,N,NB,N*N,N*NB}(translations, basissites)
    end
end
function offset_unitcell(uc::UnitCell{Sym}, offset) where Sym
    if offset isa AbstractVector{<:Number}
        @check_size offset dims(uc)
        return UnitCell{Sym}(uc.translations, uc.basissites .+ offset)
    elseif offset === :origin
        return uc
    elseif offset === :center
        basis_com = vec(sum(uc.basissites, dims=2) / length(uc))
        return UnitCell{Sym}(uc.translations, uc.basissites .- basis_com)
    elseif offset === :centeralign
        basis_com = vec(sum(uc.basissites, dims=2) / length(uc))
        uc_com = vec(sum(uc.translations, dims=2) / dims(uc))
        return UnitCell{Sym}(uc.translations, uc.basissites .- basis_com .+ uc_com)
    else
        throw(ArgumentError("invalid `offset` keyword argument"))
    end
end

"""
    UnitCell(translations[, basis; offset])

Constructs a Bravais lattice unit cell with given translation vectors and locations of basis sites.

## Arguments:
- `translations`: an `AbstractMatrix` of size `N×N` representing the translation vectors of the lattice.
- `basis`: an `AbstractMatrix` of size `N×NB` representing the locations of basis sites.
    If not provided, the lattice basis will consist of one site located in the bottom-left corner of the unit cell.

## Keyword arguments:
- `offset`: a keyword argument that specifies how to shift the lattice basis.
    Possible values:
    - `:origin`: no shift (default).
    - `:center`: shift the lattice so that the center of the basis is at the origin of the unit cell.
    - `:centeralign`: shift the lattice so that the center of the basis is at the center of the unit cell.
    - Also accepts an `AbstractVector` of size `N` to shift the lattice by a custom vector.
"""
function UnitCell(translations::AbstractMatrix,
        basissites::AbstractMatrix=zeros(size(translations, 1), 1);
        offset=:origin)
    uc_without_offset = UnitCell{:GenericBravaisLattice}(translations, basissites)
    return offset_unitcell(uc_without_offset, offset)
end

basvector(uc::UnitCell, i::Int) = uc.basissites[:, i]
unitvector(uc::UnitCell, i::Int) = uc.translations[:, i]
unitvectors(uc::UnitCell) = uc.translations

dims(::UnitCell{Sym, N} where Sym) where {N} = N
Base.length(::UnitCell{Sym,N,NB} where {Sym,N}) where {NB} = NB
Base.:(==)(b1::UnitCell, b2::UnitCell) =
    b1.translations == b2.translations && b1.basissites == b2.basissites

function print_vectors(io::IO, a::AbstractMatrix)
    indent = getindent(io)
    println(io, indent, "┌      ┐ " ^ size(a, 2))
    for i in 1:size(a, 1)
        print(io, indent)
        for j in 1:size(a, 2)
            print(io, @sprintf("│%6.3f│ ", a[i, j]))
        end
        println(io)
    end
    print(io, indent, "└      ┘ " ^ size(a, 2))
end
Base.summary(io::IO, ::UnitCell{Sym,N,NB}) where {Sym,N,NB} =
    print(io, "Unit cell of a $N-dim $Sym (", fmtnum(NB, "site"), "in basis)")
function Base.show(io::IO, ::MIME"text/plain", b::UnitCell)
    indent = getindent(io)
    print(io, indent)
    requires_compact(io) && return summary(io, b)
    if get(io, :showtitle, true)
        summary(io, b)
        println(io, ":")
    end
    io = addindent(io, 2)
    println(io, indent, "  Basis site coordinates:")
    print_vectors(io, b.basissites)
    println(io, "\n", indent, "  Translation vectors:")
    print_vectors(io, b.translations)
end

struct BravaisPointer{N}
    latcoords::SVector{N,Int}
    basindex::Int
end
dims(::BravaisPointer{N}) where N = N

@inline site_coords(uc::UnitCell, lp::BravaisPointer) =
    uc.basissites[:, lp.basindex] + uc.translations * lp.latcoords
@inline site_coords(b::UnitCell{Sym,N,1} where {Sym}, lp::BravaisPointer{N}) where {N} =
    vec(b.basissites) + b.translations * lp.latcoords

Base.:(==)(lp1::BravaisPointer, lp2::BravaisPointer) =
    lp1.basindex == lp2.basindex && lp1.latcoords == lp2.latcoords
function Base.isless(lp1::BravaisPointer, lp2::BravaisPointer)
    if lp1.latcoords == lp2.latcoords
        return isless(lp1.basindex, lp2.basindex)
    else
        return isless(lp1.latcoords, lp2.latcoords)
    end
end

"""
    BravaisSite{N, B}
A site of a `BravaisLattice{N, B}` lattice.

Fields:
- `unitcell`: a `UnitCell` object representing the lattice unit cell.
- `latcoords`: a `SVector` of size `N` representing the lattice coordinates of the site.
- `basindex`: an `Int` representing the index of the site in the lattice basis.
- `coords`: a `SVector` of size `N` representing the spatial coordinates of the site.
"""
struct BravaisSite{N,UnitcellT} <: AbstractSite{N}
    unitcell::UnitcellT
    latcoords::SVector{N,Int}
    basindex::Int
    coords::SVector{N,Float64}
    BravaisSite(lp::BravaisPointer{N}, b::UnitcellT) where {N,UnitcellT<:UnitCell} =
        new{N,UnitcellT}(b, lp.latcoords, lp.basindex, site_coords(b, lp))
end
BravaisSite(::Nothing, ::UnitCell) = NoSite()
bravaispointer(site::BravaisSite) = BravaisPointer(site.latcoords, site.basindex)

Base.show(io::IO, ::MIME"text/plain", site::BravaisSite{N, <:UnitCell{Sym}}) where {N, Sym} =
    print(io, "Site of a ", N, "-dim ", Sym, " (Bravais) @ x = $(site.coords)")

struct LatticeCoord <: SiteProperty axis::Int end
function getsiteproperty(site::BravaisSite, c::LatticeCoord)
    @assert 1 ≤ c.axis ≤ dims(site)
    return getfield(site, :latcoords)[c.axis]
end
for i in 1:32
    @eval SitePropertyAlias{$(QuoteNode(Symbol("j$i")))}() = LatticeCoord($i)
end

struct BasisIndex <: SiteProperty end
getsiteproperty(site::BravaisSite, ::BasisIndex) = getfield(site, :basindex)
SitePropertyAlias{:index}() = BasisIndex()

Base.:(==)(site1::BravaisSite, site2::BravaisSite) =
    site1.latcoords == site2.latcoords && site1.basindex == site2.basindex && site1.coords == site2.coords
Base.isless(site1::BravaisSite, site2::BravaisSite) =
    isless(bravaispointer(site1), bravaispointer(site2))
