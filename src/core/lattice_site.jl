using Logging, StaticArrays

"""
    Bravais{N, NB}
`N`-dimensional infinite Bravais lattice with `NB` sites in basis.

---
    Bravais(translation_vectors[, basis])
Constructs a Bravais lattice with given translation vectors and locations of basis sites
relative to some unit cell.
The `basis` argument can be omitted, in which case the lattice basis will consist of
one site located in the bottom-left corner of the unit cell.

`translation_vectors` argument must be an `AbstractMatrix{<:Real}` of size `N×N`,
while `basis` must also be an  abstract matrix of size `N×NB`.
"""
struct Bravais{N,NB}
    translation_vectors::SMatrix{N,N,Float64}
    basis::SMatrix{N,NB,Float64}
    function Bravais(translation_vectors::AbstractMatrix{<:Real}, basis::AbstractMatrix{<:Real},
        origin::AbstractVector{<:Real}=zeros(size(basis)[1]))
        (size(translation_vectors)[1] != size(basis)[1]) &&
            error("inconsistent dimension count (got $(size(translation_vectors)[1]), $(size(basis)[1]))")
        N, NB = size(basis)
        new{N,NB}(translation_vectors, basis .- origin)
    end
end
Bravais(translation_vectors::AbstractMatrix{<:Real}) =
    Bravais(translation_vectors, zeros((size(translation_vectors)[1], 1)))

dims(@nospecialize _::Bravais{N}) where {N} = N
Base.length(::Bravais{N,NB}) where {N,NB} = NB

struct LatticePointer{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end
dims(::LatticePointer{N}) where N = N

"""
    LatticeSite{N}
A site of a `Lattice{LatticeSym, N, NB}` lattice.

Fields:
- `unit_cell`: a set of translations along all axes representing the unit cell the site is located in.
- `basis_index`: the number of site in the lattice basis.

This type is used to iterate over all sites of a `Lattice{LatticeSym, N, NB}`.
The exact location of a `LatticeSite` can be found using the `site.coords` function.
"""
struct LatticeSite{N}
    lp::LatticePointer{N}
    coords::SVector{N,Float64}
end
dims(::LatticeSite{N}) where N = N

axis_parse_error(sym::Symbol) = error("Symbol :$sym does not correspond to any of lattice axes")
axis_parse_error(i::Int) = error("Integer $i does not correspond to any of lattice axes")
function try_parse_axis_sym(sym::Symbol)
    sym === :index && return (:b, nothing)
    sym === :x && return (:c, 1)
    sym === :y && return (:c, 2)
    sym === :z && return (:c, 3)
    axis_s = string(sym)
    axis = tryparse(Int, axis_s[2:end])
    if axis isa Int && axis > 0
        axis_s[1] == 'x' && return (:c, axis)
        axis_s[1] == 'j' && return (:j, axis)
    end
    return nothing
end
function try_parse_axis_sym(i::Int)
    i ≤ 0 && nothing
    (:c, i)
end

function get_coord(site::LatticeSite, desc)
    axtype, index = desc
    if axtype == :c && index ≤ dims(site)
        return site.coords[index]
    elseif axtype == :j && index ≤ dims(site)
        return site.lp.unit_cell[index]
    elseif axtype == :b
        return site.lp.basis_index
    else
        return nothing
    end
end

function Base.getproperty(site::LatticeSite{N}, sym::Symbol) where N
    sym === :basis_index && return site.lp.basis_index
    sym === :unit_cell && return site.lp.unit_cell
    desc = try_parse_axis_sym(sym)
    desc === nothing && return getfield(site, sym)
    val = get_coord(site, desc)
    val === nothing ? getfield(site, sym) : val
end

Base.iterate(site::LatticeSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)

cartesian_index(lp::LatticePointer) = CartesianIndex(lp.basis_index, lp.unit_cell...)
cartesian_index(site::LatticeSite) = cartesian_index(site.lp)
function LatticePointer(cind::CartesianIndex)
    tup = Tuple(cind)
    return LatticePointer(SVector(Base.tail(tup)), tup[1])
end

Base.:(==)(lp1::LatticePointer, lp2::LatticePointer) =
    lp1.basis_index == lp2.basis_index && lp1.unit_cell == lp2.unit_cell

Base.:(==)(lp::LatticePointer, site::LatticeSite) =
    site.lp == lp
Base.:(==)(site::LatticeSite, lp::LatticePointer) =
    site.lp == lp

Base.:(==)(site1::LatticeSite, site2::LatticeSite) =
    site1.lp == site2.lp && site1.coords == site2.coords

Base.isless(site1::LatticeSite, site2::LatticeSite) =
    isless(cartesian_index(site1), cartesian_index(site2))
