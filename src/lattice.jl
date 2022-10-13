using RecipesBase, LinearAlgebra, Logging, StaticArrays
import Base: length, size, copy, iterate, eltype, show, ==

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
    function Bravais(translation_vectors::AbstractMatrix{<:Real}, basis::AbstractMatrix{<:Real})
        (size(translation_vectors)[1] != size(basis)[1]) &&
            error("inconsistent dimension count (got $(size(translation_vectors)[1]), $(size(basis)[1]))")
        N, NB = size(basis)
        new{N,NB}(translation_vectors, basis)
    end
end
function Bravais(translation_vectors::AbstractMatrix{<:Real})
    Bravais(translation_vectors, zeros((size(translation_vectors)[1], 1)))
end

dims(@nospecialize _::Bravais{N}) where {N} = N
length(::Bravais{N,NB}) where {N,NB} = NB

"""
    Lattice{LatticeSym, N, NB}
A finite subset of a `Brvais{N, NB}`. `LatticeSym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Lattice` object behavior known at compile-time,
which allows to introduce various optimizations or to define specific plot recipes.

---
    Lattice(sym, sz, bvs[, mask])
Constructs a finite `Lattice{sym, N, NB}` as a subset of the `bvs` Bravais lattice.
`sz` is a `NTuple{N, Int}` which represents how many times the unit cell of `bvs` was translated by each axis - these sites form a *macro cell*.
`mask`, if defined, is a `Vector{Bool}` storing information about which of the sites from the macro cell
are actually included in the lattice, and which are not.

For example, a 3×3 square lattice with its center site excluded is represented as
`Lattice(:square, (3, 3), Bravais([1 0; 0 1]), Bool[1, 1, 1, 1, 0, 1, 1, 1, 1])`

To define a new type of lattice, create an alias for `Lattice{YourSym, YourN, YourNB}`.
Refer to the docs for detailed explanation.
"""
struct Lattice{LatticeSym,N,NB}
    lattice_size::NTuple{N,Int}
    bravais::Bravais{N,NB}
    mask::Vector{Bool}
    function Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N,NB}, mask::Vector{Bool}) where {N,NB}
        length(mask) != prod(sz) * length(bvs) &&
            error("inconsistent mask length")
        new{sym,N,NB}(sz, bvs, mask)
    end
end
Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}) where {N} =
    Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}, fill(true, prod(sz) * length(bvs)))

==(l1::T, l2::T) where {T<:Lattice} =
    (size(l1) == size(l2)) && (bravais(l1) == bravais(l2))

copy(l::Lattice{LatticeSym}) where {LatticeSym} =
    Lattice(LatticeSym, size(l), bravais(l), copy(l.mask))
size(l::Lattice) = l.lattice_size
lattice_type(::Lattice{LatticeSym}) where {LatticeSym} = LatticeSym
dims(::Lattice{LatticeSym,N} where {LatticeSym}) where {N} = N
basis_length(::Lattice{LatticeSym,N,NB} where {LatticeSym,N}) where {NB} = NB
bravais(l::Lattice) = l.bravais
length(l::Lattice) = count(l.mask)

"""
    LatticeSite{N}
A site of a `Lattice{LatticeSym, N, NB}` lattice.
Fields:
- `unit_cell`: a set of translations along all axes representing the unit cell the site is located in.
- `basis_index`: the number of site in the lattice basis.

This type is used to iterate over all sites of a `Lattice{LatticeSym, N, NB}`.
The exact location of a `LatticeSite` can be found using the `coords(lattice, site)` function.
"""
struct LatticeSite{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end

LatticeSite(tup::NTuple{N,Int}) where {N} = LatticeSite{N - 1}(SVector(tup[1:end-1]), tup[end])

Base.eltype(::Lattice{LatticeSym,N}) where {LatticeSym,N} = LatticeSite{N}

==(@nospecialize(li1::LatticeSite), @nospecialize(li2::LatticeSite)) =
    li1.basis_index == li2.basis_index && li1.unit_cell == li2.unit_cell

function iterate(l::Lattice{LatticeSym,N,NB} where {LatticeSym}) where {N,NB}
    cinds = CartesianIndex{N + 1}(1):CartesianIndex(size(l)..., NB)
    cind, st = iterate(cinds)
    index = 1
    while !l.mask[index]
        nx = iterate(cinds, cind)
        nx === nothing && return nothing
        cind, st = nx
        index += 1
    end
    LatticeSite(Tuple(cind)), (l.mask, cinds, cind, index)
end

function iterate(::Lattice, state)
    mask, cinds, cind, index = state
    nx = iterate(cinds, cind)
    nx === nothing && return nothing
    cind, st = nx
    index += 1
    while !mask[index]
        nx = iterate(cinds, cind)
        nx === nothing && return nothing
        cind, st = nx
        index += 1
    end
    LatticeSite(Tuple(cind)), (mask, cinds, cind, index)
end

"""
    coords(lattice::Lattice, site::LatticeSite) -> vector
Finds the location in space of lattice site `site` on lattice `lattice`.
"""
coords(lattice::Lattice, site::LatticeSite) =
    bravais(lattice).basis[:, site.basis_index] + bravais(lattice).translation_vectors * (site.unit_cell - SVector(size(lattice)) / 2 .- 0.5)

coords(lattice::Lattice{LatticeSym,N,1} where {LatticeSym,N}, site::LatticeSite) =
    bravais(lattice).translation_vectors * (site.unit_cell - SVector(size(lattice)) / 2 .- 0.5)

"""
    radius_vector(lattice::Lattice, site1::LatticeSite, site2::LatticeSite) -> vector
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macro cell to minimize the distance between them).
"""
function radius_vector(lattice::Lattice, site1::LatticeSite, site2::LatticeSite)
    ret_vec = coords(lattice, site1) - coords(lattice, site2)
    tr_diff = (site1.unit_cell - site2.unit_cell) ./ size(lattice)
    bravais_tr_vecs = bravais(lattice).translation_vectors
    for i in eachindex(tr_diff)
        if tr_diff[i] > 0.5
            ret_vec -= bravais_tr_vecs[:, i] * size(lattice)[i]
        elseif tr_diff[i] < -0.5
            ret_vec += bravais_tr_vecs[:, i] * size(lattice)[i]
        end
    end
    return ret_vec
end

@recipe function f(l::Lattice; show_excluded_sites=false, high_contrast=false)
    if high_contrast
        markersize := 4
        markercolor := :black
        markerstrokealpha := 1
        markerstrokestyle := :solid
        markerstrokewidth := 2
        markerstrokecolor := :white
    end
    if show_excluded_sites
        l_without_mask = copy(l)
        fill!(l_without_mask.mask, true)
        opacity := l.mask .* 0.9 .+ 0.1
        l_without_mask, nothing
    else
        l, nothing
    end
end

function collect_coords(l::Lattice)
    d = dims(l)
    pts = zeros(d, length(l))
    i = 1
    for site in l
        pts[:, i] = coords(l, site)
        i += 1
    end
    pts
end

@recipe function f(l::Lattice, v)
    label --> nothing
    show_excluded --> false
    aspect_ratio := :equal
    marker_z := v
    markerstrokewidth --> 0
    pts = collect_coords(l)
    if dims(l) == 3
        X, Y, Z = eachrow(pts)
        Xr, Yr, Zr = eachrow(round.(pts, digits=3))
        seriestype := :scatter3d
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ", ", Zr, ")")
        end
        X, Y, Z
    else
        if dims(l) == 1
            X = vec(pts)
            Y = zero(X)
        else
            X, Y = eachrow(pts[1:2, :])
        end
        Xr, Yr = eachrow(round.(pts[1:2, :], digits=3))
        seriestype --> :scatter
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ")")
        end
        if plotattributes[:seriestype] == :scatter
            X, Y
        elseif plotattributes[:seriestype] == :surface
            X, Y, v
        else
            throw(ArgumentError("unsupported series type $(plotattributes[:seriestype])"))
        end
    end
end

function show(io::IO, ::MIME"text/plain", l::Lattice{LatticeSym,N}) where {N,LatticeSym}
    print(io, "$(length(l))-site ", LatticeSym)
    if N != 1
        print(io, " lattice on ", join(size(l), "×"), " base")
    elseif !all(l.mask)
        print(io, " chain on ", size(l)[1], "-cell base")
    else
        print(io, " chain")
    end
    if N > 1 && basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

"""
    sublattice(lf::Function, l::Lattice) -> Lattice
Generates a a subset of lattice `l` by applying the `lf` function to its sites.
The `lf` function must accept two positional arguments (a `LatticeSite` and a vector with its coordinates)
and return a boolean value.
"""
function sublattice(f::Function, l::Lattice{LatticeSym}) where {LatticeSym}
    new_mask = zero(l.mask)
    new_mask[l.mask] = [f(site, coords(l, site)) for site in l]
    Lattice(LatticeSym, size(l), bravais(l), new_mask)
end

function Lattice{T,N,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

const AnyDimLattice{T,NB} = Lattice{T,N,NB} where {N}

AnyDimLattice{T,NB}(sz::Vararg{Int,N}; kw...) where {T,N,NB} =
    Lattice{T,N,NB}(sz...; kw...)

function AnyDimLattice{T,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

"""
    SquareLattice{N}
Type alias for `Lattice{:square,N,1}`.

---
    SquareLattice(sz::Int...)

Constructs a square lattice of size `sz`.
"""
const SquareLattice{N} = Lattice{:square,N,1}
function SquareLattice{N}(sz::Vararg{Int,N}) where {N}
    eye = SMatrix{N,N}(I)
    Lattice(:square, sz, Bravais(eye))
end
coords(l::SquareLattice, site::LatticeSite) =
    site.unit_cell - SVector(size(l)) / 2 .- 0.5

"""
    HoneycombLattice
Type alias for `Lattice{:honeycomb,2,2}`.

---
    HoneycombLattice(xsz::Int, ysz::Int)

Constructs a honeycomb lattice with a `xsz`×`ysz` macro cell.
"""
const HoneycombLattice = Lattice{:honeycomb,2,2}
function HoneycombLattice(xsz::Int, ysz::Int)
    bvs = Bravais([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
    Lattice(:honeycomb, (xsz, ysz), bvs)
end
