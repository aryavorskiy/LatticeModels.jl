using Logging, StaticArrays, Printf

"""
    UnitCell{Sym, N, NB}
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
    - `:center`: shift the lattice so that the center of the basis is at the origin.
    - `:centeralign`: shift the lattice so that the center of the basis is at the center of the unit cell.
    - Also accepts an `AbstractVector` of size `N` to shift the lattice by a custom vector.
"""
function UnitCell(translations::AbstractMatrix,
        basissites::AbstractMatrix=zeros(size(translations, 1), 1);
        offset=:origin)
    uc_without_offset = UnitCell{:GenericBravaisLattice}(translations, basissites)
    return offset_unitcell(uc_without_offset, offset)
end

sublatvector(uc::UnitCell, i::Int) = uc.basissites[:, i]
trvectors(uc::UnitCell) = uc.translations

dims(::UnitCell{Sym, N} where Sym) where {N} = N
Base.length(::UnitCell{Sym,N,NB} where {Sym,N}) where {NB} = NB
Base.:(==)(b1::UnitCell, b2::UnitCell) =
    b1.translations == b2.translations && b1.basissites == b2.basissites

function println_vectors(io::IO, a::AbstractMatrix)
    println(io, "┌     ┐ " ^ size(a, 2))
    for i in 1:size(a, 1)
        for j in 1:size(a, 2)
            s = @sprintf("%.3f", a[i, j])
            print(io, "│$(s)│ ")
        end
        println(io)
    end
    println(io, "└     ┘ " ^ size(a, 2))
end
function Base.show(io::IO, ::MIME"text/plain", b::UnitCell{Sym,N,NB}) where {Sym,N,NB}
    print(io, "Unit cell of a ", N, "-dim ", Sym, " (Bravais) with ",
        NB, " basis site", NB == 1 ? "" : "s")
    if !get(io, :compact, false)
        if NB == 1
            println(";")
        else
            println(io, ":")
            println_vectors(io, b.basissites)
        end
        println(io, "Translation vectors:")
        println_vectors(io, b.translations)
    end
end

struct BravaisPointer{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end
dims(::BravaisPointer{N}) where N = N

@inline site_coords(uc::UnitCell, lp::BravaisPointer) =
    uc.basissites[:, lp.basis_index] + uc.translations * lp.unit_cell
@inline site_coords(b::UnitCell{Sym,N,1} where {Sym}, lp::BravaisPointer{N}) where {N} =
    vec(b.basissites) + b.translations * lp.unit_cell

Base.:(==)(lp1::BravaisPointer, lp2::BravaisPointer) =
    lp1.basis_index == lp2.basis_index && lp1.unit_cell == lp2.unit_cell
function Base.isless(lp1::BravaisPointer, lp2::BravaisPointer)
    if lp1.unit_cell == lp2.unit_cell
        return isless(lp1.basis_index, lp2.basis_index)
    else
        return isless(lp1.unit_cell, lp2.unit_cell)
    end
end

"""
    BravaisSite{N, B}
A site of a `BravaisLattice{N, B}` lattice.

Fields:
- `lp`: a `BravaisPointer` object representing the location of the site in the Bravais lattice.
- `unit_cell`: a `UnitCell` object representing the lattice unit cell.
- `coords`: a `SVector` of size `N` representing the spatial coordinates of the site.
"""
struct BravaisSite{N,UnitcellT} <: AbstractSite{N}
    lp::BravaisPointer{N}
    unitcell::UnitcellT
    coords::SVector{N,Float64}
    BravaisSite(lp::BravaisPointer{N}, b::UnitcellT) where {N,UnitcellT<:UnitCell} =
        new{N,UnitcellT}(lp, b, site_coords(b, lp))
end
BravaisSite(::Nothing, ::UnitCell) = NoSite()

Base.show(io::IO, ::MIME"text/plain", site::BravaisSite{N, <:UnitCell{Sym}}) where {N, Sym} =
    print(io, "Site of a ", N, "-dim ", Sym, " (Bravais) @ x = $(site.coords)")

struct UnitcellAxis <: SiteProperty axis::Int end
function getsiteproperty(site::BravaisSite, c::UnitcellAxis)
    @assert 1 ≤ c.axis ≤ dims(site)
    return getfield(site, :lp).unit_cell[c.axis]
end
for i in 1:32
    @eval SitePropertyAlias{$(QuoteNode(Symbol("j$i")))}() = UnitcellAxis($i)
end

struct UnitcellIndex <: SiteProperty end
getsiteproperty(site::BravaisSite, ::UnitcellIndex) = getfield(site, :lp).basis_index
SitePropertyAlias{:index}() = UnitcellIndex()

Base.:(==)(site1::BravaisSite, site2::BravaisSite) =
    site1.lp == site2.lp && site1.coords == site2.coords
Base.isless(site1::BravaisSite, site2::BravaisSite) =
    isless(site1.lp, site2.lp)
