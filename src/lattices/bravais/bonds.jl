using StaticArrays

"""
    BravaisShift{T, N}

A struct representing bonds in some direction in a lattice.

---
    BravaisShift([site_indices, ]translate_uc)

Constructs a `BravaisShift` object.

## Arguments:
- `site_indices`: A `::Int => ::Int` pair with indices of sites connected by the bond.
When not defined, the resulting `bonds` object will connect site with any basis index to
a site with the same basis index, but in another unit cell.
- `translate_uc`: The unit cell offset.

If `site_indices` are equal or undefined and `translate_uc` is zero, the bond connects
each site with itself. In this case an error will be thrown.
Note that though the dimension count for the bond is static, it is automatically compatible with higher-dimensional lattices.
"""
struct BravaisShift{LT, N} <: AbstractTranslation{LT}
    lat::LT
    site_indices::Pair{Int, Int}
    translate_uc::SVector{N, Int}
    function BravaisShift(latt::LT, site_indices::Pair{Int, Int}, tr_uc::AbstractVector) where LT<:AbstractLattice
        if any(<(1), site_indices) && site_indices != (0 => 0)
            throw(ArgumentError("Invalid site indices $site_indices: ≥1 or 0=>0 expected"))
        end
        iszero(tr_uc) && ==(site_indices...) && throw(ArgumentError("bond connects site to itself"))
        new{LT, length(tr_uc)}(latt, site_indices, tr_uc)
    end
end
BravaisShift(lat::AbstractLattice, tr_uc::AbstractVector) = BravaisShift(lat, 0=>0, tr_uc)
BravaisShift(args...; kw...) = BravaisShift(UndefinedLattice(), args...; kw...)
apply_lattice(bsh::BravaisShift{UndefinedLattice}, l::AbstractLattice) =
    BravaisShift(l, bsh.site_indices, bsh.translate_uc)
dims(::BravaisShift{UndefinedLattice, N}) where N = N

@inline has_sublatremap(bsh::BravaisShift) = bsh.site_indices != (0 => 0)

"""
    BravaisShift(site_indices)
    BravaisShift([site_indices; ]axis[, dist=1])

A convenient constructor for a `BravaisShift` object.

## Arguments:
- `site_indices`: a `::Int => ::Int` pair with indices of sites connected by the bond;
if omitted, the bond connects sites with the same sublattice index.

## Keyword arguments:
- `axis`: The hopping direction axis in terms of unit cell vectors.
- `dist`: The hopping distance in terms of
"""
function BravaisShift(lat::AbstractLattice, site_indices::Pair{Int,Int} = 0=>0; axis=0, dist=1)
    axis == 0 && return BravaisShift(site_indices, [])
    BravaisShift(lat, site_indices, one_hot(axis, axis) * dist)
end

Base.:(==)(h1::BravaisShift, h2::BravaisShift) =
    all(getfield(h1, fn) == getfield(h2, fn) for fn in fieldnames(BravaisShift))
function Base.inv(bsh::BravaisShift)
    a, b = bsh.site_indices
    BravaisShift(bsh.lat, b => a, -bsh.translate_uc)
end

function Base.show(io::IO, mime::MIME"text/plain", bsh::BravaisShift)
    println(io, "BravaisShift; Unit cell shift ", bsh.translate_uc,
        ", Sublattice remapping: ",
        has_sublatremap(bsh) ? "none" : bsh.site_indices)
    if !(bsh.lat isa UndefinedLattice)
        print(io, "\n on")
        show(io, mime, bsh.lat)
    end
end

@inline function _destination_bp(bsh::BravaisShift{LT, N} where LT, lp::BravaisPointer{M}) where {N, M}
    if has_sublatremap(bsh)
        bsh.site_indices[1] != lp.basis_index && return nothing
        new_basindex = bsh.site_indices[2]
    else
        new_basindex = lp.basis_index
    end
    N > M && any(!=(0), @view bsh.translate_uc[M+1:N]) && return nothing
    return BravaisPointer(add_assuming_zeros(lp.unit_cell, bsh.translate_uc), new_basindex)
end
@inline destination(bs::BravaisShift, site::BravaisSite) =
    BravaisSite(_destination_bp(bs, site.lp), site.bravais)
