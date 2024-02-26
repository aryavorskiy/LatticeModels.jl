using StaticArrays

"""
    BravaisTranslation{T, N}

A struct representing bonds in some direction in a lattice.

---
Note that though the dimension count for the bond is static, it is automatically compatible with higher-dimensional lattices.
"""
struct BravaisTranslation{LT, N} <: AbstractTranslation{LT}
    lat::LT
    site_indices::Pair{Int, Int}
    translate_uc::SVector{N, Int}
    function BravaisTranslation(latt::LT, site_indices::Pair{Int, Int}, tr_uc::AbstractVector) where LT<:AbstractLattice
        if any(<(1), site_indices) && site_indices != (0=>0)
            throw(ArgumentError("Invalid site indices $site_indices: ≥1 expected"))
        end
        iszero(tr_uc) && ==(site_indices...) && throw(ArgumentError("bond connects site to itself"))
        new{LT, length(tr_uc)}(latt, site_indices, SVector{length(tr_uc)}(tr_uc))
    end
end
BravaisTranslation(lat::AbstractLattice, tr_uc::AbstractVector) = BravaisTranslation(lat, 0=>0, tr_uc)
BravaisTranslation(args...; kw...) = BravaisTranslation(UndefinedLattice(), args...; kw...)
adapt_bonds(bsh::BravaisTranslation, l::AbstractLattice) =
    BravaisTranslation(l, bsh.site_indices, bsh.translate_uc)
dims(::BravaisTranslation{UndefinedLattice, N}) where N = N

struct Bravais end
Base.getindex(::Type{Bravais}, I::Int...) = BravaisTranslation(UndefinedLattice(), 0=>0, SVector(I))

"""
    BravaisTranslation([site_indices, ]translate_uc)
    BravaisTranslation(site_indices)
    BravaisTranslation([site_indices; ]axis[, dist=1])

A convenient constructor for a `BravaisTranslation` object.

## Arguments:
- `site_indices`: a `::Int => ::Int` pair with indices of sites connected by the bond;
if omitted, the bond connects sites with the same sublattice index.
- `translate_uc`: The unit cell offset.

## Keyword arguments:
- `axis`: The hopping direction axis in terms of unit cell vectors.
- `dist`: The hopping distance in terms of unit cell vectors.

If `site_indices` are equal or undefined and `translate_uc` is zero, the translation is
considered to be a translation of all sites to themselves. An error will be thrown in this case.
"""
function BravaisTranslation(lat::AbstractLattice, site_indices::Pair{Int,Int} = 0=>0; axis=0, dist=1)
    axis == 0 && return BravaisTranslation(site_indices, [])
    BravaisTranslation(lat, site_indices, one_hot(axis, axis) * dist)
end

Base.:(==)(h1::BravaisTranslation, h2::BravaisTranslation) =
    all(getfield(h1, fn) == getfield(h2, fn) for fn in fieldnames(BravaisTranslation))
function Base.inv(bsh::BravaisTranslation)
    a, b = bsh.site_indices
    BravaisTranslation(bsh.lat, b => a, -bsh.translate_uc)
end

function Base.summary(io::IO, bsh::BravaisTranslation)
    if bsh.site_indices == (0=>0)
        print(io, "Bravais", bsh.translate_uc)
    else
        print(io, "$(bsh.site_indices), $(bsh.translate_uc)")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", bsh::BravaisTranslation)
    indent = getindent(io)
    get(io, :showtitle, true) && !requires_compact(io) &&
        print(io, indent, "BravaisTranslation:\n")
    print(io, indent)
    summary(io, bsh)
    if !(bsh.lat isa UndefinedLattice) && !requires_compact(io)
        print(io, "\n", indent, " on ")
        show(io, mime, bsh.lat)
    end
end

@inline function _destination_bp(bsh::BravaisTranslation{LT, N} where LT, lp::BravaisPointer{M}) where {N, M}
    if bsh.site_indices == (0=>0)
        new_basindex = lp.basindex
    else
        bsh.site_indices[1] != lp.basindex && return nothing
        new_basindex = bsh.site_indices[2]
    end
    N > M && any(!=(0), @view bsh.translate_uc[M+1:N]) && return nothing
    return BravaisPointer(add_assuming_zeros(lp.latcoords, bsh.translate_uc), new_basindex)
end
@inline destination(bs::BravaisTranslation, site::BravaisSite) =
    BravaisSite(_destination_bp(bs, bravaispointer(site)), site.unitcell)

struct BravaisSiteMapping{LT, TupleT} <: DirectedBonds{LT}
    lat::LT
    translations::TupleT
    function BravaisSiteMapping(latt::LT, translations::Vararg{BravaisTranslation{UndefinedLattice}}) where LT<:AbstractLattice
        new{LT, typeof(translations)}(latt, translations)
    end
end
BravaisSiteMapping(trs::BravaisTranslation{UndefinedLattice}...) =
    BravaisSiteMapping(UndefinedLattice(), trs...)
adapt_bonds(bts::BravaisSiteMapping, l::AbstractLattice) =
    BravaisSiteMapping(l, bts.translations...)

Base.inv(bsm::BravaisSiteMapping) = BravaisSiteMapping(bsm.lat, inv.(bsm.translations)...)
@inline destinations(bsm::BravaisSiteMapping, site::BravaisSite) =
    (destination(tr, site) for tr in bsm.translations)

merge_lats(l1::AbstractLattice, l2::AbstractLattice) = throw(IncompatibleLattices(l1, l2))
merge_lats(l::AbstractLattice, ::UndefinedLattice) = l
merge_lats(::UndefinedLattice, l::AbstractLattice) = l

function Base.union(tr1::BravaisTranslation, tr2::BravaisTranslation)
    u_tr1 = adapt_bonds(tr1, UndefinedLattice())
    u_tr2 = adapt_bonds(tr2, UndefinedLattice())
    l = merge_lats(tr1.lat, tr2.lat)
    if u_tr1 == u_tr2
        return BravaisSiteMapping(l, u_tr1)
    else
        return BravaisSiteMapping(l, u_tr1, u_tr2)
    end
end
function Base.union(trs::BravaisSiteMapping, tr::BravaisTranslation)
    u_tr = adapt_bonds(tr, UndefinedLattice())
    l = merge_lats(trs.lat, tr.lat)
    if u_tr in trs.translations
        return trs
    else
        return BravaisSiteMapping(l, trs.translations..., u_tr)
    end
end
Base.union(tr::BravaisTranslation, trs::BravaisSiteMapping) = union(trs, tr)
Base.union(trs::BravaisSiteMapping, trs2::BravaisSiteMapping) =
    foldl(union, trs2.translations, init=trs)

function Base.summary(io::IO, trs::BravaisSiteMapping)
    print(io, "BravaisSiteMapping with $(length(trs.translations)) translations")
end
function Base.show(io::IO, ::MIME"text/plain", trs::BravaisSiteMapping)
    if get(io, :showtitle, true)
        print(io, getindent(io))
        summary(io, trs)
        print(io, ": ")
        requires_compact(io) || println(io)
        io = addindent(io)
    end
    if requires_compact(io)
        get(io, :showtitle, true) || print(io, getindent(io))
        return print(io, "($(fmtnum(trs.translations, "translation")) not shown)")
    else
        for i in 1:length(trs.translations)
            print(io, getindent(io))
            if get(io, :maxlines, typemax(Int)) ≤ i < length(trs.translations)
                print(io, " ⋮ ")
                omitted = length(trs.translations) - i + 1
                get(io, :showtitle, true) || print(io, "($omitted more omitted)")
                break
            end
            summary(io, trs.translations[i])
            i < length(trs.translations) && println(io)
        end
    end
    if !(trs.lat isa UndefinedLattice)
        print(io, "\n on ")
        summary(io, trs.lat)
    end
end

function adapt_bonds(tr::Translation{UndefinedLattice}, l::OnSites{BravaisLattice{N}}) where N
    M = length(tr.R)
    if M > N && !all(==(0), @view tr.R[N+1:M])
        throw(ArgumentError("Cannot adapt a $M-dim translation to a $N-dim lattice"))
    end
    shifts = BravaisTranslation{UndefinedLattice, N}[]
    for i in 1:baslength(l)
        for j in 1:baslength(l)
            R′ = add_assuming_zeros(basvector(l, i) - basvector(l, j), tr.R)
            J = unitvectors(l) \ R′
            if all(<(√eps()), abs.(rem.(J, 1, RoundNearest)))
                i == j && return BravaisTranslation(J)
                # J is an integer vector
                push!(shifts, BravaisTranslation(i=>j, J))
                break
            end
        end
    end
    return BravaisSiteMapping(l, shifts...)
end
Translation(l::OnSites{BravaisLattice}, R::AbstractVector) = adapt_bonds(Translation(R), l)

function detect_nnhops(uc::UnitCell{Sym, N,NB} where Sym, depth=2, limit=5) where {N,NB}
    lens = Float64[]
    TranslationT = BravaisTranslation{UndefinedLattice, N}
    translation_lists = Vector{TranslationT}[]
    for c in cartesian_indices(depth, Val(N))
        C = SVector{N}(Tuple(c))
        for i in 1:length(uc)
            for j in i:length(uc)
                if i == j
                    nzc = findfirst(!=(0), C)
                    nzc === nothing && continue # skip self-hops
                    C[nzc] < 1 && continue # avoid double-counting
                end
                R = unitvectors(uc) * C - basvector(uc, i) + basvector(uc, j)
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
    return DefaultNNBonds(Tuple(lens[1:3]), Tuple(trs[1:3]))
end
