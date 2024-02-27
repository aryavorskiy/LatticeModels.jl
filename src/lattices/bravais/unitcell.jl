using Logging, StaticArrays, Printf

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
struct UnitCell{N,NU,NB,NND,NNB}
    translations::SMatrix{N,NU,Float64,NND}
    basissites::SMatrix{N,NB,Float64,NNB}
    function UnitCell(translations::AbstractMatrix,
            basissites::AbstractMatrix=zeros(size(translations, 1), 1);
            offset=:origin)
        N, NU = size(translations)
        NB = size(basissites, 2)
        @check_size basissites (N, NB)
        u = new{N,NU,NB,N*NU,N*NB}(translations, basissites)
        return offset_unitcell(u, offset)
    end
end
function offset_unitcell(uc::UnitCell, offset)
    if offset isa AbstractVector{<:Number}
        @check_size offset dims(uc)
        return UnitCell(uc.translations, uc.basissites .+ offset)
    elseif offset === :origin
        return uc
    elseif offset === :center
        basis_com = vec(sum(uc.basissites, dims=2) / length(uc))
        return UnitCell(uc.translations, uc.basissites .- basis_com)
    elseif offset === :centeralign
        basis_com = vec(sum(uc.basissites, dims=2) / length(uc))
        uc_com = vec(sum(uc.translations, dims=2) / dims(uc))
        return UnitCell(uc.translations, uc.basissites .- basis_com .+ uc_com)
    else
        throw(ArgumentError("invalid `offset` keyword argument"))
    end
end

basvector(uc::UnitCell, i::Int) = uc.basissites[:, i]
unitvector(uc::UnitCell, i::Int) = uc.translations[:, i]
unitvectors(uc::UnitCell) = uc.translations

dims(::UnitCell{N}) where {N} = N
ldims(::UnitCell{N,NU}) where {N,NU} = NU
Base.length(::UnitCell{N,NU,NB}) where {N,NU,NB} = NB
Base.:(==)(b1::UnitCell, b2::UnitCell) =
    b1.translations == b2.translations && b1.basissites == b2.basissites
lattransform(ltr::LatticeTransform, uc::UnitCell) = UnitCell(uc.translations |> ltr, uc.basissites |> ltr)

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
Base.summary(io::IO, ::UnitCell{N,ND,NB}) where {N,ND,NB} =
    print(io, "$NB-site Unit cell of a $ND-dim Bravais lattice in $(N)D space")
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

struct BravaisPointer{NU}
    latcoords::SVector{NU,Int}
    basindex::Int
end
ldims(::BravaisPointer{NU}) where NU = NU

@inline site_coords(uc::UnitCell, lp::BravaisPointer) =
    uc.basissites[:, lp.basindex] + uc.translations * lp.latcoords
@inline site_coords(b::UnitCell{N,NU,1}, lp::BravaisPointer{NU}) where {N,NU} =
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

lattransform(::LatticeTransform, lp::BravaisPointer) = lp

"""
    BravaisSite{N,NU,B}
A site of a `BravaisLattice{N,NU,B}` lattice.

Fields:
- `unitcell`: a `UnitCell` object representing the lattice unit cell.
- `latcoords`: a `SVector` of size `N` representing the lattice coordinates of the site.
- `basindex`: an `Int` representing the index of the site in the lattice basis.
- `coords`: a `SVector` of size `N` representing the spatial coordinates of the site.
"""
struct BravaisSite{N,NU,UnitcellT} <: AbstractSite{N}
    unitcell::UnitcellT
    latcoords::SVector{NU,Int}
    basindex::Int
    coords::SVector{N,Float64}
    BravaisSite(lp::BravaisPointer{NU}, b::UnitcellT) where {N,NU,UnitcellT<:UnitCell{N,NU}} =
        new{N,NU,UnitcellT}(b, lp.latcoords, lp.basindex, site_coords(b, lp))
end
BravaisSite(::Nothing, ::UnitCell) = NoSite()
unitcell(site::BravaisSite) = site.unitcell
ldims(::BravaisSite{N,NU,UnitcellT}) where {N,NU,UnitcellT} = NU
bravaispointer(site::BravaisSite) = BravaisPointer(site.latcoords, site.basindex)
lattransform(ltr::LatticeTransform, site::BravaisSite) =
    BravaisSite(bravaispointer(site) |> ltr, site.unitcell |> ltr)

struct LatticeCoord <: SiteProperty axis::Int end
function getsiteproperty(site::BravaisSite, c::LatticeCoord)
    1 ≤ c.axis ≤ ldims(site) ||
        throw(ArgumentError("invalid lattice axis $(c.axis) for $(ldims(site))-dimensional Bravais site"))
    return getfield(site, :latcoords)[c.axis]
end
for i in 1:32
    @eval SitePropertyAlias{$(QuoteNode(Symbol("j$i")))}() = LatticeCoord($i)
end

lattransform(proj::Project{LatticeCoord}, uc::UnitCell{N,NU}) where {N, NU} =
    UnitCell(uc.translations[:, proj.projection.axis .!= 1:NU], uc.basissites)
lattransform(proj::Project{LatticeCoord}, lp::BravaisPointer{NU}) where NU =
    BravaisPointer(SVector{NU-1}(lp.latcoords[proj.projection.axis .!= 1:NU]), lp.basindex)

struct BasisIndex <: SiteProperty end
getsiteproperty(site::BravaisSite, ::BasisIndex) = getfield(site, :basindex)
SitePropertyAlias{:index}() = BasisIndex()

lattransform(::Project{BasisIndex}, uc::UnitCell{N}) where {N} =
    UnitCell(uc.translations, zero(SVector{N}))
lattransform(::Project{BasisIndex}, lp::BravaisPointer) =
    BravaisPointer(lp.latcoords, 1)

Base.:(==)(site1::BravaisSite, site2::BravaisSite) =
    site1.latcoords == site2.latcoords && site1.basindex == site2.basindex && site1.coords == site2.coords
Base.isless(site1::BravaisSite, site2::BravaisSite) =
    isless(bravaispointer(site1), bravaispointer(site2))

sitekey(site::BravaisSite) = site.latcoords[1]
secondarykey(site::BravaisSite{N,2} where N) = site.basindex + site.latcoords[2] * length(site.unitcell)
