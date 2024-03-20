"""
    NearestNeighbor{N}

A bonds type that connects sites that are nearest neighbors of order `N` on some lattice.
"""
struct NearestNeighbor{N} <: AbstractBonds{UndefinedLattice}
    function NearestNeighbor(::Val{N}) where {N}
        new{N}()
    end
end
NearestNeighbor(N::Int=1) = NearestNeighbor(Val(N))

struct DefaultNNBonds{M, TupleT}
    dists::NTuple{M, Float64}
    nnbonds::TupleT
    function DefaultNNBonds(dists::NTuple{M, Float64}, nnbonds::Tuple) where {M}
        length(dists) == length(nnbonds) ||
            throw(ArgumentError("The number of distances and the number of hops must be the same."))
        all(isfinite, dists) || all(isnan, dists) ||
            throw(ArgumentError("All distances must be finite"))
        nnbonds_nolat = Tuple(adapt_bonds(b, UndefinedLattice()) for b in nnbonds)
        new{M, typeof(nnbonds_nolat)}(dists, nnbonds_nolat)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", dnn::DefaultNNBonds)
    indent = getindent(io)
    print(io, indent, "Nearest neighbor hoppings: ")
    isempty(dnn.dists) && return print(io, "none")
    io = addindent(io, 2, :showtitle => false)
    for i in 1:length(dnn.dists)
        if isfinite(dnn.dists[i])
            println(io, "\n", indent, @sprintf("%9.5f", dnn.dists[i]), " =>")
        else
            println(io, "\n", indent, "  #$i =>")
        end
        show(io, mime, dnn.nnbonds[i])
    end
end

Base.getindex(dnn::DefaultNNBonds, i::Int) = dnn.nnbonds[i]
Base.length(dnn::DefaultNNBonds) = length(dnn.nnbonds)

"""
    getnnbonds(lat)

Returns the nearest neighbor bonds of the lattice `lat`.
"""
getnnbonds(l::AbstractLattice) = getmeta(l, :nnbonds, DefaultNNBonds((), ()))
setnnbonds(l::AbstractLattice, dnn::DefaultNNBonds) = setmeta(l, :nnbonds, dnn)
_to_nnb(b::AbstractBonds) = NaN => b
_to_nnb(p::Pair{Float64,<:AbstractBonds}) = p
_to_nnb(b) =
    throw(ArgumentError("Invalid nnbond type: $(typeof(b)); expected a bonds type or a distance-bonds pair"))

"""
    setnnbonds(lat, args...; overwrite=false)

Adds the nearest neighbor bonds `args` to the lattice `lat`. If `overwrite` is `true`, the default
nearest neighbor bonds are replaced by `args`. Otherwise, the new bonds are merged with the default.

Each `args` can be a bonds type or a distance-bonds pair.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(3, 3);

julia> lat2 = setnnbonds(lat, SiteDistance(0..1), SiteDistance(1..2));

julia> lat2.nnbonds
Nearest neighbor hoppings:
  #1 =>
    SiteDistance(0 .. 1)
  #2 =>
    SiteDistance(1 .. 2)
```
"""
function setnnbonds(l::AbstractLattice, args::AbstractBonds...)
    new_nn = DefaultNNBonds(Tuple(NaN for p in args), args)
    return setnnbonds(l, new_nn)
end

adapt_bonds(::NearestNeighbor, ::AbstractLattice) = NoBonds()
function adapt_bonds(b::NearestNeighbor{N}, l::LatticeWithMetadata) where {N}
    default_nnhops = getnnbonds(l)
    if default_nnhops === nothing || N > length(default_nnhops)
        return adapt_bonds(adapt_bonds(b, l.lat), l)
    else
        return adapt_bonds(default_nnhops[N], l)
    end
end
NearestNeighbor(l::LatticeWithMetadata, N=1) = adapt_bonds(NearestNeighbor(N), l)

function detect_nnhops(uc::UnitCell{Sym, N,NB} where Sym, depth=2, limit=3) where {N,NB}
    lens = Float64[]
    TranslationT = BravaisTranslation{UndefinedLattice, N}
    translation_lists = Vector{TranslationT}[]
    for c in cartesian_indices(depth, Val(N))
        C = SVector{N}(Tuple(c))
        U = unitvectors(uc) * C
        for i in 1:length(uc)
            for j in i:length(uc)
                if i == j
                    i > 2 && continue # this was already done for i == j == 1
                    nzc = findfirst(!=(0), C)
                    nzc === nothing && continue # skip self-hops
                    C[nzc] < 1 && continue # avoid double-counting
                end
                R = U - basvector(uc, i) + basvector(uc, j)
                r = norm(R)
                k = findfirst(l -> lens[l] > r || lens[l] ≈ r, eachindex(lens))
                if k === nothing || !(lens[k] ≈ r)
                    k === nothing && (k = length(lens) + 1)
                    k > limit && continue
                    insert!(lens, k, r)
                    insert!(translation_lists, k, TranslationT[])
                end

                # if i == j, the bond is valid for all basis indices
                translation = i == j ? BravaisTranslation(C) : BravaisTranslation(i=>j, C)
                if !(translation in translation_lists[k])
                    push!(translation_lists[k], translation)
                end
            end
        end
    end
    return lens, [BravaisSiteMapping(UndefinedLattice(), trs...) for trs in translation_lists]
end

function adapt_bonds(::NearestNeighbor{N}, l::BravaisLattice) where N
    adapt_bonds(detect_nnhops(l.unitcell, 1 + ceil(Int, √N), N)[2][N], l)
end

"""
    NearestNeighbor(lat[, N=1])

Returns the nearest neighbor bonds of order `N` for the lattice `lat`.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = HoneycombLattice(5, 5);

julia> NearestNeighbor(lat)
BravaisSiteMapping with 3 translations:
  1 => 2, [0, -1]
  1 => 2, [-1, 0]
  1 => 2, [0, 0]
 on 50-site 2-dim Bravais lattice in 2D space (2-site basis)

julia> lat = SquareLattice(3, 3, 3, 3);

julia> NearestNeighbor(lat, 4)
BravaisSiteMapping with 12 translations:
  Bravais[1, -1, -1, -1]
  Bravais[1, 1, -1, -1]
  Bravais[1, -1, 1, -1]
  Bravais[1, 1, 1, -1]
  Bravais[2, 0, 0, 0]
  Bravais[0, 2, 0, 0]
  Bravais[0, 0, 2, 0]
  Bravais[1, -1, -1, 1]
  Bravais[1, 1, -1, 1]
   ⋮
 on 81-site 4-dim Bravais lattice in 4D space
```
"""
NearestNeighbor(l::BravaisLattice, N=1) = adapt_bonds(NearestNeighbor(N), l)

function getnnbonds(l::BravaisLattice)
    lens, trs = detect_nnhops(l.unitcell)
    l = min(3, length(lens))
    return DefaultNNBonds(Tuple(lens[1:l]), Tuple(trs[1:l]))
end
