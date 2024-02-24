abstract type Boundary{TranslationT} end
function phasefactor end

function nshifts(site::AbstractSite, tr::AbstractTranslation, n::Int)
    n == 0 && site
    n < 0 && (tr = -tr)
    for _ in 1:abs(n)
        site -= tr
    end
    site
end
function nshifts_phase(site::AbstractSite, b::Boundary, n::Int)
    n == 0 && return 1., site
    factor = 1.
    for _ in 1:abs(n)
        if n > 0
            site -= b.translation
            factor *= phasefactor(b, site)
        else
            factor /= phasefactor(b, site)
            site += b.translation
        end
    end
    factor, site
end

to_boundary(any) = throw(ArgumentError("Could not convert `$any` to a `Boundary`"))
to_boundary(::Nothing) = nothing
to_boundary(b::Boundary) = b

"""
    TwistedBoundary <: Boundary

A boundary condition with a phase twist. A `PeriodicBoundary` is a special case of `TwistedBoundary` with zero twist.

---
    TwistedBoundary(translation, Θ)
    PeriodicBoundary(translation)

Construct a `TwistedBoundary` with a given translation and twist angle.

## Arguments
- `translation`: The translation vector of the boundary representad as `AbstractTranslation`.
    If an array is passed, it is converted to `Translation` automatically.
- `Θ`: The twist angle in radians.
"""
struct TwistedBoundary{TranslationT} <: Boundary{TranslationT}
    translation::TranslationT
    Θ::Float64
    TwistedBoundary(tr::TranslationT, Θ::Real) where TranslationT<:AbstractTranslation =
        new{TranslationT}(tr, Float64(Θ))
end
TwistedBoundary(tr::AbstractArray, Θ::Real) = TwistedBoundary(Translation(tr), Θ)
PeriodicBoundary(tr) = TwistedBoundary(tr, 0)
adapt_boundary(b::TwistedBoundary, l::AbstractLattice) = TwistedBoundary(adapt_bonds(b.translation, l), b.Θ)

to_boundary(p::Pair{<:Any, Bool}) =
    p[2] ? PeriodicBoundary(p[1]) : nothing
to_boundary(p::Pair{<:Any, <:Real}) =
    TwistedBoundary(p[1], p[2])
function Base.show(io::IO, mime::MIME"text/plain", bc::TwistedBoundary)
    show(io, mime, bc.translation)
    print(io, " → ")
    if bc.Θ % 2π ≈ 0
        print(io, "periodic")
    else
        print(io, "twist θ = ", trunc(bc.Θ, digits=2))
    end
end

phasefactor(b::TwistedBoundary, ::AbstractSite) = exp(im * b.Θ)

@doc raw"""
    FunctionBoundary <: Boundary

A boundary condition with a function that returns the phase factor for a given site.
The boundary condition is encoded in form $ψ(x + R) = f(x)ψ(x)$, where $f(x)$ is the
function and $R$ is the translation vector.

---
    FunctionBoundary(f, translation)

Construct a `FunctionBoundary` with a given function and translation.

## Arguments
- `f`: The function that returns the phase factor for a given site.
- `translation`: The translation vector of the boundary representad as `AbstractTranslation`.
    If an array is passed, it is converted to `Translation` automatically.
"""
struct FunctionBoundary{TranslationT<:AbstractTranslation, F<:Function} <: Boundary{TranslationT}
    condition::F
    translation::TranslationT
end
FunctionBoundary(f::F, tr::AbstractArray) where F<:Function =
    FunctionBoundary(f, Translation(tr))
adapt_boundary(b::FunctionBoundary, l::AbstractLattice) =
    FunctionBoundary(b.condition, adapt_bonds(b.translation, l))
to_boundary(p::Pair{<:Any, <:Function}) = FunctionBoundary(p[2], p[1])
function Base.show(io::IO, mime::MIME"text/plain", bc::FunctionBoundary)
    show(io, mime, bc.translation)
    print(io, " → function '", bc.condition, "'")
end

phasefactor(b::FunctionBoundary, site::AbstractSite) = b.condition(site)

"""
    BoundaryConditions

A collection of boundary conditions for a lattice.

## Fields
- `bcs`: A tuple of boundary conditions.
- `depth`: The upper limit of the depth of the boundary conditions (used for routing).
"""
struct BoundaryConditions{CondsTuple}
    bcs::CondsTuple
    depth::Int
    function BoundaryConditions(bcs::CondsTuple; depth=1) where CondsTuple<:Tuple{Vararg{Boundary}}
        new{CondsTuple}(bcs, depth)
    end
end

_skipnothing(args::Tuple) = _skipnothing((), args)
_skipnothing(pre_args::Tuple, post_args::Tuple) =
    first(post_args) === nothing ?
    _skipnothing(pre_args, Base.tail(post_args)) :
    _skipnothing((pre_args..., first(post_args)), Base.tail(post_args))
_skipnothing(pre_args::Tuple, ::Tuple{}) = pre_args

BoundaryConditions(args...; kw...) =
    BoundaryConditions(_skipnothing(to_boundary.(args)); kw...)

parse_translation(l::AbstractLattice, b::Boundary) = adapt_boundary(b, l)
parse_translation(l::AbstractLattice, pair::Pair) = parse_translation(l, pair[1]) => pair[2]
parse_translation(l::AbstractLattice, tr::AbstractTranslation) = adapt_bonds(tr, l)
parse_translation(l::AbstractLattice, sym::Symbol) = defaulttranslations(l)[sym]
parse_translation(l::AbstractLattice, vec::AbstractVector) = adapt_bonds(Translation(vec), l)
parse_translation(::AbstractLattice, any) = throw(ArgumentError("Could not interpret `$any` as a Translation"))

function parse_boundaries(l::AbstractLattice, arg)
    if arg isa Boundary
        return BoundaryConditions(arg)
    elseif arg isa BoundaryConditions
        return arg
    elseif arg isa Tuple
        return BoundaryConditions(parse_translation.(Ref(l), arg)...)
    else
        return BoundaryConditions(parse_translation(l, arg))
    end
end

adapt_boundaries(bcs::BoundaryConditions, l::AbstractLattice) =
    BoundaryConditions(map(b -> adapt_boundary(b, l), bcs.bcs); depth=bcs.depth)
getboundaries(l::AbstractLattice) = adapt_boundaries(getparam(l, :boundaries, BoundaryConditions()), l)
setboundaries(l::AbstractLattice, bcs::BoundaryConditions) =
    setparam(l, :boundaries, adapt_boundaries(bcs, UndefinedLattice()))
setboundaries(l::AbstractLattice, bcs) = setboundaries(l, parse_boundaries(l, bcs))

Base.getindex(bcs::BoundaryConditions, i::Int) = bcs.bcs[i]

function Base.show(io::IO, mime::MIME"text/plain", bcs::BoundaryConditions)
    indent = getindent(io)
    print(io, indent, "Boundary conditions")
    length(bcs.bcs) == 0 && return print(io, ": none")
    if requires_compact(io)
        print(io, ": (", length(bcs.bcs),
        " not shown; depth = ", bcs.depth, ")")
    else
        print(io, " (depth = ", bcs.depth, "):")
        io = addindent(io, :compact => true)
        for i in 1:length(bcs.bcs)
            println(io)
            show(io, mime, bcs[i])
        end
    end
end

cartesian_indices(b::BoundaryConditions{<:Tuple{Vararg{Any,M}}}) where M =
    cartesian_indices(b.depth, Val{M}())
cartesian_indices(l::AbstractLattice) = cartesian_indices(getboundaries(l))

function route(bcs::BoundaryConditions, l::AbstractLattice, site::AbstractSite)
    rs = resolve_site_default(l, site)
    if rs !== nothing
        tup = Tuple(zero(SVector{length(bcs.bcs), Int}))
        return tup, rs
    end
    for cind in cartesian_indices(bcs)
        tup = Tuple(cind)
        all(==(0), tup) && continue
        new_site = site
        new_site === NoSite() && continue
        for i in eachindex(tup)
            new_site = nshifts(new_site, bcs[i].translation, tup[i])
            new_site === NoSite() && break
        end
        rs = resolve_site_default(l, new_site)
        rs !== nothing && return tup, rs
    end
    return nothing
end
function route(::BoundaryConditions{Tuple{}}, l::AbstractLattice, site::AbstractSite)
    rs = resolve_site_default(l, site)
    rs === nothing ? nothing : ((), rs)
end

function findfactor(site::AbstractSite, tup, bounds::Tuple{Vararg{Boundary}})
    factor = 1.0 + 0.0im
    for i in 1:length(tup)
        new_factor, site = nshifts_phase(site, bounds[i], tup[i])
        factor *= new_factor
    end
    return factor
end
@inline function findfactor(::AbstractSite, tup, bounds::Tuple{Vararg{TwistedBoundary}})
    return exp(im * sum(bounds[i].Θ * tup[i] for i in 1:length(tup)))
end

function resolve_site(l::LatticeWithParams, site::AbstractSite)
    bcs = getboundaries(l)
    r = route(bcs, l, site)
    r === nothing && return nothing
    tup, rs = r
    all(==(0), tup) && return rs
    factor = findfactor(site, tup, bcs.bcs)
    return ResolvedSite(rs.site, site, rs.index, factor * rs.factor)
end

"""
    translate_to_nearest(lat, site1, site2)

Translate `site2` to its equivalent nearest to `site1` in the lattice `lat`, taking the
boundary conditions into account.
"""
function translate_to_nearest(l::AbstractLattice, site1::AbstractSite, site2::AbstractSite)
    min_site = site2
    min_dist = norm(site1.coords - site2.coords)
    for cind in cartesian_indices(l)
        tup = Tuple(cind)
        new_site = site2
        for i in eachindex(tup)
            new_site = nshifts(new_site, getboundaries(l)[i].translation, tup[i])
        end
        if norm(new_site.coords - site1.coords) < min_dist
            min_dist = norm(new_site.coords - site1.coords)
            min_site = new_site
        end
    end
    return min_site
end
site_distance(l::LatticeWithParams, site1::AbstractSite, site2::AbstractSite) =
    norm(translate_to_nearest(l, site1, site2).coords - site1.coords)

struct DefaultTranslations{NamedTupleT}
    translations::NamedTupleT
end
function DefaultTranslations(translations::Vararg{Pair{Symbol, <:AbstractTranslation}})
    ntup = NamedTuple(translations)
    DefaultTranslations(ntup)
end
function DefaultTranslations(dt::DefaultTranslations, translations::Pair{Symbol, <:AbstractTranslation}...)
    ntup = merge(dt.translations, translations)
    DefaultTranslations(ntup)
end
function Base.getindex(dt::DefaultTranslations, sym::Symbol)
    sym in keys(dt.translations) ||
        throw(ArgumentError("No translation for symbol :$sym. One of the following is available: $(keys(dt.translations))"))
    return dt.translations[sym]
end

function Base.show(io::IO, mime::MIME"text/plain", dt::DefaultTranslations)
    print(io, "Default translations: ")
    isempty(dt.translations) && return print(io, "none")
    io = IOContext(io, :compact => true)
    for (sym, tr) in pairs(dt.translations)
        print(io, "\n  :", sym, " → ")
        show(io, mime, tr)
    end
end

function addtranslations(l::AbstractLattice, translations::Pair{Symbol, <:AbstractTranslation}...; overwrite=false)
    tr = getparam(l, :defaulttranslations, DefaultTranslations())
    if overwrite
        setparam(l, :defaulttranslations, DefaultTranslations(translations...))
    else
        setparam(l, :defaulttranslations, DefaultTranslations(tr, translations...))
    end
end
addtranslations(l::AbstractLattice, translations::Tuple; kw...) =
    addtranslations(l, translations...; kw...)

function defaulttranslations(l::AbstractLattice)
    getparam(l, :defaulttranslations, DefaultTranslations())
end
