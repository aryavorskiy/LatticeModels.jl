"""
    NearestNeighbor{N}

A bonds type that connects sites that are nearest neighbors of order `N` on some lattice.
"""
struct NearestNeighbor{N} <: AbstractBonds{UndefinedLattice}
    function NearestNeighbor(::Val{N}) where {N}
        new{N}()
    end
end
NearestNeighbor(N::Int) = NearestNeighbor(Val(N))

struct DefaultNNBonds{M, TupleT}
    dists::NTuple{M, Float64}
    nnbonds::TupleT
    function DefaultNNBonds(dists::NTuple{M, Float64}, nnbonds::Tuple) where {M}
        length(dists) == length(nnbonds) ||
            throw(ArgumentError("The number of distances and the number of hops must be the same."))
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
        println(io, "\n", indent, @sprintf("%9.5f", dnn.dists[i]), " =>")
        show(io, mime, dnn.nnbonds[i])
    end
end

Base.getindex(dnn::DefaultNNBonds, i::Int) = dnn.nnbonds[i]
Base.length(dnn::DefaultNNBonds) = length(dnn.nnbonds)

"""
    getnnbonds(lat)

Returns the nearest neighbor bonds of the lattice `lat`.
"""
getnnbonds(l::AbstractLattice) = getparam(l, :nnbonds, DefaultNNBonds((), ()))
setnnbonds(l::AbstractLattice, dnn::DefaultNNBonds) = setparam(l, :nnbonds, dnn)

function adapt_bonds(b::NearestNeighbor{N}, l::LatticeWithParams) where {N}
    default_nnhops = getnnbonds(l)
    if default_nnhops === nothing || N > length(default_nnhops)
        return adapt_bonds(adapt_bonds(b, l.lat), l)
    else
        return adapt_bonds(default_nnhops[N], l)
    end
end
NearestNeighbor(l::LatticeWithParams, N) = adapt_bonds(NearestNeighbor(N), l)

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
NearestNeighbor(l::BravaisLattice, N) = adapt_bonds(NearestNeighbor(N), l)

function getnnbonds(l::BravaisLattice)
    lens, trs = detect_nnhops(l.unitcell)
    l = min(3, length(lens))
    return DefaultNNBonds(Tuple(lens[1:l]), Tuple(trs[1:l]))
end
