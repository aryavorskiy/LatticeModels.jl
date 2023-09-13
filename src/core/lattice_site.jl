using Logging, StaticArrays

"""
    Bravais{Sym, N, NB}
`N`-dimensional infinite Bravais lattice with `NB` sites in basis.
`Sym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Bravais` object behavior known at compile-time,
which allows to introduce various optimizations or define specific plot recipes.

---
    Bravais(translation_vectors[, basis, origin])
Constructs a Bravais lattice with given translation vectors and locations of basis sites
relative to some unit cell.
The `basis` argument can be omitted, in which case the lattice basis will consist of
one site located in the bottom-left corner of the unit cell.

`translation_vectors` argument must be an `AbstractMatrix{<:Real}` of size `N×N`,
while `basis` must also be an  abstract matrix of size `N×NB`.
"""
struct Bravais{Sym,N,NB}
    translation_vectors::SMatrix{N,N,Float64}
    basis::SMatrix{N,NB,Float64}
    function Bravais{Sym}(translation_vectors::AbstractMatrix{<:Real},
        basis::AbstractMatrix{<:Real}=zeros((size(translation_vectors, 1), 1)),
        origin::AbstractVector{<:Real}=zeros(size(basis)[1])) where Sym
        (size(translation_vectors)[1] != size(basis)[1]) &&
            error("inconsistent dimension count (got $(size(translation_vectors)[1]), $(size(basis)[1]))")
        N, NB = size(basis)
        new{Sym,N,NB}(translation_vectors, basis .+ origin)
    end
end

dims(@nospecialize _::Bravais{Sym, N} where Sym) where {N} = N
Base.:(==)(b1::BT, b2::BT) where BT<:Bravais =
    b1.translation_vectors == b2.translation_vectors && b1.basis == b2.basis
Base.length(::Bravais{Sym,N,NB} where {Sym,N}) where {NB} = NB

struct LatticePointer{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end
dims(::LatticePointer{N}) where N = N

site_coords(b::Bravais, lp::LatticePointer) =
    b.basis[:, lp.basis_index] + b.translation_vectors * lp.unit_cell
site_coords(b::Bravais{Sym,N,1} where {Sym}, lp::LatticePointer{N}) where {N} =
    vec(b.basis) + b.translation_vectors * lp.unit_cell

function Base.isless(site1::LatticePointer, site2::LatticePointer)
    if site1.unit_cell == site2.unit_cell
        return isless(site1.basis_index, site2.basis_index)
    else
        return isless(site1.unit_cell, site2.unit_cell)
    end
end

"""
    LatticeSite{N}
A site of a `Lattice{LatticeSym, N, NB}` lattice.

Fields:
- `unit_cell`: a set of translations along all axes representing the unit cell the site is located in.
- `basis_index`: the number of site in the lattice basis.

This type is used to iterate over all sites of a `Lattice{LatticeSym, N, NB}`.
The exact location of a `LatticeSite` can be found using the `site.coords` function.
"""
struct LatticeSite{N, B}
    lp::LatticePointer{N}
    bravais::B
    coords::SVector{N,Float64}
    LatticeSite(lp::LatticePointer{N}, b::B) where {N, B<:Bravais} =
        new{N, B}(lp, b, site_coords(b, lp))
end
dims(::LatticeSite{N}) where N = N

const AxisSpec = Tuple{Symbol, Nullable{Int}}
function try_parse_axis_sym(sym::Symbol)
    sym === :index && return (:b, nothing)
    sym === :basis_index && return (:b, nothing)
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
try_parse_axis_sym(i::Int) = i > 0 ? (:c, i) : nothing
function parse_axis_sym(sym)
    desc = try_parse_axis_sym(sym)
    desc === nothing && throw(ArgumentError("Invalid axis specifier '$sym'"))
    return desc
end

function get_coord(site::LatticeSite, desc::AxisSpec)
    axtype, index = desc
    if axtype == :c && index ≤ dims(site)
        return site.coords[index]
    elseif axtype == :j && index ≤ dims(site)
        return site.lp.unit_cell[index]
    elseif axtype == :b && index === nothing
        return site.lp.basis_index
    elseif index > dims(site)
        throw(DimensionMismatch("$N-dimensional site does not have axis #$index"))
    else
        error("Invalid coord specifier $desc")
    end
end

function Base.getproperty(site::LatticeSite{N}, sym::Symbol) where N
    sym == :unit_cell && return site.lp.unit_cell
    desc = try_parse_axis_sym(sym)
    if desc === nothing
        return getfield(site, sym)
    else
        return get_coord(site, desc)
    end
end

Base.iterate(site::LatticeSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)

Base.:(==)(lp1::LatticePointer, lp2::LatticePointer) =
    lp1.basis_index == lp2.basis_index && lp1.unit_cell == lp2.unit_cell
Base.:(==)(site1::LatticeSite, site2::LatticeSite) =
    site1.lp == site2.lp && site1.coords == site2.coords
