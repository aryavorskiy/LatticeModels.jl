using RecipesBase, LinearAlgebra, Logging, StaticArrays
import Base: length, size, copy, iterate, eltype, show, ==

struct Bravais{N,NB}
    translation_vectors::SMatrix{N,N,Float64}
    basis::SMatrix{N,NB,Float64}
    function Bravais(translation_vectors::AbstractMatrix{<:Real}, basis::AbstractMatrix{<:Real})
        (size(translation_vectors)[1] != size(basis)[1]) && error("inconsistent dimension count")
        N, NB = size(basis)
        new{N,NB}(translation_vectors, basis)
    end
end
function Bravais(translation_vectors::AbstractMatrix{<:Real})
    Bravais(translation_vectors, zeros(Float64, (size(translation_vectors)[1], 1)))
end

dims(@nospecialize _::Bravais{N}) where {N} = N
length(::Bravais{N,NB}) where {N,NB} = NB

struct Lattice{LatticeType,N}
    lattice_size::NTuple{N,Int}
    bravais::Bravais{N}
    mask::Vector{Bool}
    function Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}, mask::Vector{Bool}) where {N}
        length(mask) != prod(sz) * length(bvs) &&
            error("inconsistent mask length")
        new{sym,N}(sz, bvs, mask)
    end
end
Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}) where {N} =
    Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}, fill(true, prod(sz) * length(bvs)))

==(@nospecialize(l1::Lattice), @nospecialize(l2::Lattice)) =
    (l1.lattice_size == l2.lattice_size) && (l1.bravais == l2.bravais)

size(@nospecialize l::Lattice) = l.lattice_size
dims(@nospecialize _::Lattice{LatticeType,N}) where {LatticeType,N} = N
bravais(@nospecialize l::Lattice) = l.bravais

length(@nospecialize l::Lattice) = count(l.mask)
struct LatticeIndex{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
end

Base.eltype(::Lattice{LatticeType,N}) where {LatticeType,N} = LatticeIndex{N}

==(@nospecialize(li1::LatticeIndex), @nospecialize(li2::LatticeIndex)) =
    li1.basis_index == li2.basis_index && li1.unit_cell == li2.unit_cell
_first(l::Lattice) = LatticeIndex(@SVector(fill(1, dims(l))), 1)

function _next(cinds::CartesianIndices, bvs_len, site::LatticeIndex)
    if site.basis_index == bvs_len
        n = iterate(cinds, CartesianIndex(Tuple(site.unit_cell)))
        if n === nothing
            nothing
        else
            LatticeIndex(SVector(Tuple(n[1])), 1)
        end
    else
        LatticeIndex(site.unit_cell, site.basis_index + 1)
    end
end

function iterate(@nospecialize l::Lattice{LatticeType,N}) where {LatticeType,N}
    site = _first(l)
    cinds = CartesianIndex{N}(1):CartesianIndex(size(l))
    bl = length(l.bravais)
    index = 1
    while !l.mask[index]
        site = _next(cinds, bl, site)
        site === nothing && return nothing
        index += 1
    end
    (site, (site, cinds, bl, index))
end

function iterate(@nospecialize(l::Lattice), state)
    site, cinds, bl, index = state
    site = _next(cinds, bl, site)
    site === nothing && return nothing
    index += 1
    while !l.mask[index]
        site = _next(cinds, bl, site)
        site === nothing && return nothing
        index += 1
    end
    (site, (site, cinds, bl, index))
end

coords(l::Lattice, site::LatticeIndex) =
    bravais(l).basis[:, site.basis_index] + bravais(l).translation_vectors * (site.unit_cell - SVector(size(l)) / 2 .- 0.5)

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

function Lattice{T}(f::Function, sz::Int...) where {T}
    l = Lattice{T}(sz...)
    sublattice(f, l)
end

function Lattice{:square}(sz::Vararg{Int,N}) where {N}
    eye = SMatrix{N,N}(I)
    Lattice(:square, sz, Bravais(eye))
end

const SquareLattice = Lattice{:square}

function Lattice{:honeycomb}(xsz::Int, ysz::Int)
    bvs = Bravais([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
    Lattice(:honeycomb, (xsz, ysz), bvs)
end

const HoneycombLattice = Lattice{:honeycomb}
