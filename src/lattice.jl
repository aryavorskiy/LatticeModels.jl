using RecipesBase, LinearAlgebra, Logging, StaticArrays
import Base: length, size, copy, iterate, eltype, show, ==

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

struct Lattice{LatticeType,N,NB}
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

==(@nospecialize(l1::Lattice), @nospecialize(l2::Lattice)) =
    (l1.lattice_size == l2.lattice_size) && (l1.bravais == l2.bravais)

size(l::Lattice) = l.lattice_size
dims(::Lattice{LatticeType,N}) where {LatticeType,N} = N
bravais(l::Lattice) = l.bravais

length(l::Lattice) = count(l.mask)
struct LatticeIndex{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end

LatticeIndex(tup::NTuple{N,Int}) where {N} = LatticeIndex{N - 1}(SVector(tup[1:end-1]), tup[end])

Base.eltype(::Lattice{LatticeType,N}) where {LatticeType,N} = LatticeIndex{N}

==(@nospecialize(li1::LatticeIndex), @nospecialize(li2::LatticeIndex)) =
    li1.basis_index == li2.basis_index && li1.unit_cell == li2.unit_cell

function iterate(l::Lattice{LT,N,NB} where {LT}) where {N,NB}
    site = LatticeIndex(@SVector(fill(1, N)), 1)
    cinds = CartesianIndex{N + 1}(1):CartesianIndex(l.lattice_size..., NB)
    cind, st = iterate(cinds)
    index = 1
    while !l.mask[index]
        nx = iterate(cinds, cind)
        nx === nothing && return nothing
        cind, st = nx
        index += 1
    end
    (LatticeIndex(Tuple(cind)), (l.mask, cinds, cind, index))
end

function iterate(@nospecialize(_::Lattice), state)
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
    (LatticeIndex(Tuple(cind)), (mask, cinds, cind, index))
end

coords(l::Lattice, site::LatticeIndex) =
    bravais(l).basis[:, site.basis_index] + bravais(l).translation_vectors * (site.unit_cell - SVector(size(l)) / 2 .- 0.5)

coords(l::Lattice{LT,N,1} where {LT,N}, site::LatticeIndex) =
    bravais(l).translation_vectors * (site.unit_cell - SVector(size(l)) / 2 .- 0.5)

function radius_vector(l::Lattice, site1::LatticeIndex, site2::LatticeIndex)
    ret_vec = coords(l, site1) - coords(l, site2)
    tr_diff = (site1.unit_cell - site2.unit_cell) ./ size(l)
    bravais_tr_vecs = bravais(l).translation_vectors
    for i in eachindex(tr_diff)
        if tr_diff[i] > 0.5
            ret_vec -= bravais_tr_vecs[:, i] * size(l)[i]
        elseif tr_diff[i] < -0.5
            ret_vec += bravais_tr_vecs[:, i] * size(l)[i]
        end
    end
    return ret_vec
end

@recipe function f(l::Lattice, v=nothing)
    label --> nothing
    show_excluded --> false
    # TODO support show_excluded=true
    aspect_ratio := :equal
    marker_z := v
    markerstrokewidth := 0
    d = dims(l)
    pts = zeros(Float64, d, length(l))
    i = 1
    for site in l
        crd = coords(l, site)
        pts[:, i] = crd
        i += 1
    end
    if dims(l) == 1
        pts = vcat(pts, zeros(Float64, 1, length(l)))
    end
    if dims(l) == 3
        X, Y, Z = eachrow(pts)
        Xr, Yr, Zr = eachrow(round.(pts, digits=3))
        seriestype := :scatter3d
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ", ", Zr, ")")
        end
        X, Y, Z
    else
        X, Y = eachrow(pts[1:2, :])
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

function show(io::IO, ::MIME"text/plain", l::Lattice{LatticeType,N}) where {N,LatticeType}
    print(io, "$(length(l))-site ", LatticeType)
    if N != 1
        print(io, " lattice on ", join(size(l), "×"), " base")
    elseif !all(l.mask)
        print(io, " chain on ", size(l)[1], "-cell base")
    else
        print(io, " chain")
    end
    if N > 1 && length(l.bravais) > 1
        print(io, " (", length(l.bravais), "-site basis)")
    end
end

function _propagate_lattice_args(f, l::Lattice)
    if hasmethod(f, NTuple{dims(l),Number})
        (_l::Lattice, site::LatticeIndex) -> f(coords(_l, site)...)
    elseif hasmethod(f, Tuple{LatticeIndex})
        (::Lattice, site::LatticeIndex) -> f(site)
    elseif hasmethod(f, Tuple{LatticeIndex,Vararg{Number,dims(l)}})
        (_l::Lattice, site::LatticeIndex) -> f(site, coords(_l, site)...)
    else
        throw(ArgumentError("failed to propagate args: unsupported lambda type"))
    end
end

_always_true_on_lattice(::Lattice, ::LatticeIndex) = true

function sublattice(f::Function, l::Lattice{LatticeType}) where {LatticeType}
    lf = _propagate_lattice_args(f, l)
    Lattice(LatticeType, size(l), bravais(l),
        [lf(l, site) for site in l])
end

function Lattice{T,N,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

const SizableLattice{T,NB} = Lattice{T,N,NB} where {N}

SizableLattice{T,NB}(sz::Vararg{Int,N}; kw...) where {T,N,NB} =
    Lattice{T,N,NB}(sz...; kw...)

function SizableLattice{T,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

const SquareLattice{N} = Lattice{:square,N,1}
function SquareLattice{N}(sz::Vararg{Int,N}) where {N}
    eye = SMatrix{N,N}(I)
    Lattice(:square, sz, Bravais(eye))
end
coords(l::SquareLattice, site::LatticeIndex) =
    site.unit_cell - SVector(size(l)) / 2 .- 0.5

const HoneycombLattice = Lattice{:honeycomb,2,2}
function HoneycombLattice(xsz::Int, ysz::Int)
    bvs = Bravais([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
    Lattice(:honeycomb, (xsz, ysz), bvs)
end
