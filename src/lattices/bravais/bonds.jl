using StaticArrays

"""
    SiteOffset{T, N}

A struct representing bonds in some direction in a lattice.

---
    SiteOffset([site_indices, ]translate_uc)

Constructs a `SiteOffset` object.

## Arguments:
- `site_indices`: A `::Int => ::Int` pair with indices of sites connected by the bond.
When not defined, the resulting `bonds` object will connect site with any basis index to
a site with the same basis index, but in another unit cell.
- `translate_uc`: The unit cell offset.

If `site_indices` are equal or undefined and `translate_uc` is zero, the bond connects
each site with itself. In this case an error will be thrown.
Note that though the dimension count for the bond is static, it is automatically compatible with higher-dimensional lattices.
"""
struct SiteOffset{T, N}
    site_indices::T
    translate_uc::SVector{N, Int}
    function SiteOffset(site_indices::Pair{Int, Int}, tr_uc::AbstractVector)
        any(<(1), site_indices) && throw(ArgumentError("Positive site indices expected"))
        iszero(tr_uc) && ==(site_indices...) && throw(ArgumentError("bond connects site to itself"))
        new{Pair{Int, Int}, length(tr_uc)}(site_indices, tr_uc)
    end
    function SiteOffset(tr_uc::AbstractVector)
        iszero(tr_uc) && throw(ArgumentError("bond connects site to itself"))
        new{Nothing, length(tr_uc)}(nothing, tr_uc)
    end
end
SiteOffset(::Nothing, tr_uc) = SiteOffset(tr_uc)

"""
    SiteOffset(site_indices)
    SiteOffset([site_indices; ]axis[, dist=1])

A convenient constructor for a `SiteOffset` object.

## Arguments:
- `site_indices`: a `::Int => ::Int` pair with indices of sites connected by the bond;
if omitted, the bond connects sites with the same sublattice index.

## Keyword arguments:
- `axis`: The hopping direction axis in terms of unit cell vectors.
- `dist`: The hopping distance in terms of
"""
function SiteOffset(site_indices::Nullable{Pair{Int,Int}}=nothing; axis=0, dist=1)
    axis == 0 && return SiteOffset(site_indices, [])
    SiteOffset(site_indices, one_hot(axis, axis) * dist)
end
const Bonds{N} = SiteOffset{N}

Base.:(==)(h1::SiteOffset, h2::SiteOffset) =
    all(getproperty(h1, fn) == getproperty(h2, fn) for fn in fieldnames(SiteOffset))

function Base.show(io::IO, ::MIME"text/plain", hop::SiteOffset{<:Pair})
    println(io, "SiteOffset connecting site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.translate_uc)")
end
function Base.show(io::IO, ::MIME"text/plain", hop::SiteOffset{<:Nothing})
    println(io, "SiteOffset connecting sites translated by $(hop.translate_uc)")
end
dims(::SiteOffset{T, N} where T) where N = N

@inline function Base.:(+)(lp::BravaisPointer{N}, bs::SiteOffset{<:Pair}) where N
    bs.site_indices[1] != lp.basis_index && return nothing
    any(bs.translate_uc[N+1:end] .!= 0) && return nothing
    return BravaisPointer(add_assuming_zeros(lp.unit_cell, bs.translate_uc), bs.site_indices[2])
end
@inline Base.:(+)(lp::BravaisPointer, bs::SiteOffset{<:Nothing}) =
    return BravaisPointer(add_assuming_zeros(lp.unit_cell, bs.translate_uc), lp.basis_index)

@inline Base.:(+)(site::BravaisSite, bs) = BravaisSite(site.lp + bs, site.bravais)

function adjacency_matrix(l::AbstractLattice{<:BravaisSite}, bss::SiteOffset...)
    matrix = zeros(Bool, length(l), length(l))
    for bs in bss
        for (i, site) in enumerate(l)
            offset_site = site + bs
            _, new_site = shift_site(l, offset_site)
            j = site_index(l, new_site)
            j === nothing && continue
            matrix[i, j] = matrix[j, i] = true
        end
    end
    return AdjacencyMatrix(l, matrix)
end
