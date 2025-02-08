import IntervalSets: Interval, leftendpoint, rightendpoint, width, mean

# N-dim unit sphere volume
ndvol(N) = N ≤ 1 ? N + 1 : ndvol(N - 2) * 2pi / N
# Scale the unit vectors to make the unit cell hold a unit sphere
function inscribe_sphere(uc::UnitCell)
    _orth(a, b) = a - b * dot(a, b) / dot(b, b)
    mv = MVector{dims(uc), Float64}(undef)
    for i in 1:dims(uc)
        v = unitvector(uc, i)
        for j in 1:dims(uc)
            j == i && continue
            w = unitvector(uc, j)
            for k in 1:j - 1
                k == i && continue
                w = _orth(w, unitvector(uc, k))
            end
            v = _orth(v, w)
        end
        mv[i] = 1 / norm(v)
    end
    return Tuple(mv)
end

abstract type AbstractShape{N} end
function circumscribed_sphere end
function inshape end
dims(::AbstractShape{N}) where N = N
function bounding_region(uc::UnitCell, shape::AbstractShape)
    r, c = circumscribed_sphere(uc, shape)
    mds = inscribe_sphere(uc)
    c_b = ceil.(Int, unitvectors(uc) \ c)
    return Tuple(-ceil(Int, r * mds[i]) + c_b[i] - 1:ceil(Int, r * mds[i]) + c_b[i]
        for i in eachindex(mds))
end

function scalefactor(uc::UnitCell, shapes::AbstractShape...; sites::Int)
    n_estimate = sum(volume, shapes) * length(uc) / det(unitvectors(uc))
    n_estimate == 0 && throw(ArgumentError("Cannot rescale shapes with zero volume"))
    return (sites / n_estimate) ^ (1 / dims(uc))
end
scalefactor(::Type{LT}, shapes::AbstractShape...; sites::Int) where LT<:BravaisLatticeType =
    scalefactor(construct_unitcell(LT), shapes...; sites=sites)

Base.:(*)(shape::AbstractShape, c::Real) = scale(shape, c)

struct NotShape{N, S} <: AbstractShape{N}
    shape::S
    NotShape(s::AbstractShape) = new{dims(s), typeof(s)}(s)
end
inshape(f::NotShape, site) = !inshape(f.shape, site)
volume(f::NotShape) = -volume(f.shape)
scale(f::NotShape, c) = NotShape(scale(f.shape, c))
Base.:(!)(s::AbstractShape) = NotShape(s)
Base.:(!)(n::NotShape) = n.shape


"""
    shaperadius(unitcell, shape, sites)
    shaperadius(lat, shape[, sites])

Calculate the radius of a shape such that it contains appriximately `sites` sites.

## Arguments
- `unitcell`: The `UnitCell` of the lattice. Might also be a lattice type.
- `lat`: The lattice. It is considered that the lattice was constructed in the same shape.
- `shape`: The shape to calculate the radius for.
- `sites`: The number of sites the shape should contain.
"""
shaperadius(uc::UnitCell, shape::AbstractShape, sites::Int) =
    shape.radius * scalefactor(uc, shape; sites=sites)
shaperadius(LT::Type{<:BravaisLatticeType}, shape::AbstractShape, sites::Int) =
    shaperadius(construct_unitcell(LT), shape, sites)
shaperadius(lat::MaybeWithMetadata{BravaisLattice}, shape::AbstractShape, sites=length(lat)) =
    shaperadius(lat.unitcell, shape, sites)

function insertordered!(arr, arr2)
    isempty(arr2) && return
    if isempty(arr)
        append!(arr, arr2)
        return
    end

    i = 1
    for j in eachindex(arr2)
        while i ≤ length(arr) && arr[i] < arr2[j]
            i += 1
        end
        if i > length(arr)
            append!(arr, @view arr2[j:end])
            break
        end
        arr[i] == arr2[j] && continue
        insert!(arr, i, arr2[j])
    end
end

"""
    fillshapes(uc, shapes...[; sites, scale, kw...])

Create a lattice sample with geometry defined by the given shapes. The lattice is filled with
sites that are inside the shapes.

## Arguments
- `uc`: The `UnitCell` of the lattice. Might also be a lattice type.
- `shapes`: The shapes to fill the lattice with.

## Keyword Arguments
- `sites`: If given, an attepmt will be made to fill the lattice with the given number of sites.
    The scaling will be approximate and relying on assumptions that the shapes do not overlap.
- `scale`: The scaling factor for the shapes. If `sites` is given, the scaling factor will be
    calculated automatically.
All other keyword arguments are passed to the lattice constructor. See [`span_unitcells`](@ref)
for more information.
"""
function fillshapes(uc::UnitCell{Sym,N} where Sym, shapes::AbstractShape...;
        sites::Nullable{Int}=nothing, scale::Real=1, offset=:origin, rotate=nothing, removedangling=2, kw...) where N
    bps = BravaisPointer{N}[]
    bps_buffer = BravaisPointer{N}[]
    if sites !== nothing
        scale != 1 && @warn "Ignoring scale factor when `sites` is given"
        scale = scalefactor(uc, shapes...; sites=sites)
        new_shapes = LatticeModels.scale.(shapes, scale)
    else
        new_shapes = shapes
    end
    new_unitcell = transform_unitcell(uc, rotate=rotate, offset=offset)
    for shape in new_shapes
        dims(shape) != N &&
            throw(ArgumentError("$(dims(shape))-dim $(typeof(shape)) incompatible with $N-dim lattice"))
        if shape isa NotShape
            filter!(bps) do bp
                inshape(shape, BravaisSite(bp, new_unitcell))
            end
        else
            empty!(bps_buffer)
            add_bravaispointers!(site -> inshape(shape, site), bps_buffer, uc, bounding_region(uc, shape))
            insertordered!(bps, bps_buffer)
        end
    end
    b = BravaisLattice(new_unitcell, bps)
    fb = finalize_lattice(b; kw...)
    removedangling!(fb, removedangling)
    return addtranslations(fb, overwrite=true)
end
fillshapes(LT::Type{<:BravaisLatticeType}, shapes::AbstractShape...; kw...) =
    settype(fillshapes(construct_unitcell(LT), shapes...; kw...), LT)

(::Type{LT})(shapes::AbstractShape...; kw...) where LT<:BravaisLatticeType =
    settype(fillshapes(construct_unitcell(LT), shapes...; kw...), LT)

function addshapes!(l::MaybeWithMetadata{BravaisLattice{N}}, shapes::AbstractShape...) where N
    bps = l.pointers
    uc = l.unitcell
    for shape in shapes
        shape isa NotShape && throw(ArgumentError("Shape negation is not supported in `addshapes!`"))
        dims(shape) != N &&
            throw(ArgumentError("$(dims(shape))-dim $(typeof(shape)) incompatible with $N-dim lattice"))
        add_bravaispointers!(site -> inshape(shape, site), bps, uc, bounding_region(uc, shape))
    end
    return l
end

function deleteshapes!(l::AbstractLattice, shapes::AbstractShape...)
    filter!(site -> !any(inshape.(shapes, Ref(site))), l)
end

"""
    BallND{N}([radius, center])

Construct a `N`-dimensional ball with a given radius and center. Note the aliases: `Circle`
and `Ball` are `BallND{2}` and `BallND{3}` respectively.

## Arguments
- `radius`: The radius of the ball.
- `center`: The center of the ball.
"""
struct BallND{N} <: AbstractShape{N}
    radius::Float64
    center::SVector{N, Float64}
    function BallND{N}(radius::Real, center::AbstractVector{<:Real}=zero(SVector{N})) where N
        @check_size center N
        new{length(center)}(radius, center)
    end
end
BallND{N}(center::AbstractVector{<:Real}=zero(SVector{N})) where N = BallND{N}(1, center)
function circumscribed_sphere(uc::UnitCell, f::BallND)
    return f.radius, add_assuming_zeros(zero(SVector{dims(uc)}), f.center)
end
scale(f::BallND{N}, c) where N = BallND{N}(c * f.radius, c * f.center)
volume(f::BallND{N}) where N = ndvol(dims(f)) * f.radius^N
inshape(f::BallND, site::AbstractSite) = sum(abs2, site.coords - f.center) ≤ f.radius^2
topath2d(f::BallND{2}) =
    [Tuple(f.center .+ f.radius * SA[cos(θ), sin(θ)]) for θ in range(0, 2pi, length=100)]

const Circle = BallND{2}
const Ball = BallND{3}

"""
    Polygon{N}([radius, center])
    Polygon{N}([center; h])

Construct a regular `N`-sided polygon with a given (circumscribed) radius and center. Note
the aliases: `Triangle`, `Square`, and `Hexagon` are `Polygon{3}`, `Polygon{4}`, and
`Polygon{6}` respectively.

## Arguments
- `radius`: The (circumscribed) radius of the polygon.
- `center`: The center of the polygon.

## Keyword Arguments
- `h`: The distance from the center to the vertices. If given, the `radius` is calculated as
  `h / cos(pi / N)`.
"""
struct Polygon{N} <: AbstractShape{2}
    radius::Float64
    center::SVector{2, Float64}
    function Polygon{N}(radius::Real, center::AbstractVector{<:Real}=zero(SVector{2})) where N
        @check_size center 2
        new{N}(radius, center)
    end
end
Polygon{N}(center::AbstractVector{<:Real}=zero(SVector{2}); h=cos(pi / N)) where N =
    Polygon{N}(h / cos(pi / N) + 1e-8, center)
function circumscribed_sphere(::UnitCell, f::Polygon)
    return f.radius, f.center
end
scale(f::Polygon{N}, c) where N = Polygon{N}(c * f.radius, c * f.center)
volume(f::Polygon{N}) where N = N / 2 * sin(2pi / N) * f.radius^2
@generated function inshape(f::Polygon{N}, site::AbstractSite) where N
    exprs = Expr[:(y * $(cos(2pi * i / N)) + x * $(sin(2pi * i / N)) > -h) for i in 1:N]
    and_expr = reduce((x, y) -> :($x && $y), exprs)
    quote
        h = f.radius * $(cos(pi / N))
        x, y = site.coords - f.center
        return $and_expr
    end
end
topath2d(f::Polygon{N}) where N =
    [Tuple(f.center .+ f.radius * SA[cos(2pi * i / N), sin(2pi * i / N)]) for i in 0:N]

const Triangle = Polygon{3}
const Square = Polygon{4}
const Hexagon = Polygon{6}

"""
    SiteAt(coords)

Represents a single site at the given coordinates.
"""
struct SiteAt{N} <: AbstractShape{N}
    coords::SVector{N, Int}
    SiteAt(center::AbstractVector{<:Int}) = new{length(center)}(center)
end
circumscribed_sphere(::UnitCell, f::SiteAt) = 1, f.coords
scale(::SiteAt, _) = throw(ArgumentError("SiteAt does not support scaling"))
volume(::SiteAt) = 0
inshape(f::SiteAt, site::AbstractSite) = isapprox(site.coords, f.coords, atol=1e-10)
topath2d(f::SiteAt) = [Tuple(f.coords)]

"""
    Box(intervals...)

Construct a box with given horizontal and vertical intervals. Usage:
`Box(1 .. 3, 2 .. 4)`.

## Arguments
- `intervals`: The intervals for each dimension.
"""
struct Box{N, T} <: AbstractShape{N}
    intervals::NTuple{N, ClosedInterval{T}}
end
function Box(ints::Vararg{Interval, N}) where N
    T = promote_type(eltype.(ints)...)
    Box{N, T}(ClosedInterval.(ints))
end
function circumscribed_sphere(::UnitCell, f)
    return sqrt(sum(abs2 ∘ width, f.intervals)) / 2,
        SVector(Tuple(mean(int) for int in f.intervals))
end
scale(int::Interval{L, R}, c) where {L, R} =
    Interval{L, R}(c * leftendpoint(int), c * rightendpoint(int))
scale(f::Box, c) = Box(Tuple(scale(int, c) for int in f.intervals))
volume(f::Box) = prod(width, f.intervals)
inshape(f::Box{N}, site::AbstractSite) where N = all(i -> site.coords[i] in f.intervals[i], SOneTo{N}())
function topath2d(f::Box{2})
    x1, x2 = endpoints(f.intervals[1])
    y1, y2 = endpoints(f.intervals[2])
    return [(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)]
end

"""
    Path(start, stop)

Construct a path from `start` to `stop`.

## Arguments
- `start`: The start of the path.
- `stop`: The end of the path.
"""
struct Path{N} <: AbstractShape{N}
    start::SVector{N, Float64}
    stop::SVector{N, Float64}
    Path(start::AbstractVector{<:Number}, stop::AbstractVector{<:Number}) =
        new{length(start)}(start, stop)
end
function circumscribed_sphere(::UnitCell, f::Path)
    return norm(f.start - f.stop)/2, (f.start + f.stop) / 2
end
scale(f::Path, c) = Path(c * f.start, c * f.stop)
volume(::Path) = 0
function inshape(f::Path, site::BravaisSite)
    sb = 0.
    eb = 1.
    j = unitvectors(site.unitcell) \ (f.stop - f.start)
    k1 = site.latcoords .- unitvectors(site.unitcell) \ f.start
    k2 = site.latcoords  .+ 1 .- vec(unitvectors(site.unitcell) \ f.start)
    for i in 1:dims(f)
        p1, p2 = extrema((k1[i], k2[i]))
        if j[i] == 0
            if p1 ≤ 0 < p2
                continue
            else
                return false
            end
        end
        if j[i] < 0
            p1, p2 = p2, p1
        end
        sb = max(sb, p1 / j[i])
        eb = min(eb, p2 / j[i])
    end
    orth = k1 * dot(j, j) - j * dot(k1, j)
    iszero(orth) && return true
    e = orth[findfirst(!=(0), orth)]
    return sb < eb || (sb ≈ eb && e > 0)
end
topath2d(f::Path{2}) = [Tuple(f.start), Tuple(f.stop)]

@inline _iter_neighbor_indices(lat::AbstractLattice, nns::AbstractBonds, i) =
    Iterators.map(adjacentsites(nns, lat[i])) do site
        rs = resolve_site(lat, site)
        return rs === nothing ? 0 : rs.index
    end

function _add_counts!(counts, bs::AbstractBonds)
    for b in bs
        b === nothing && continue
        s1, s2 = b
        s1 === nothing && continue
        s2 === nothing && continue
        counts[s1.index] += 1
        counts[s2.index] += 1
    end
end
function _add_counts!(counts, bsm::BravaisSiteMapping)
    l = lattice(bsm)
    for t in bsm.translations
        _add_counts!(counts, adapt_bonds(t, l))
    end
end

function mergesorted!(arr, arr2)
    i = 1
    for j in eachindex(arr2)
        i = searchsortedfirst(arr, arr2[j], i, length(arr), Base.Order.Forward)
        if i > length(arr)
            append!(arr, @view arr2[j:end])
            break
        end
        insert!(arr, i, arr2[j])
        i += 1
    end
end

function _removedangling!(rlat, maxdepth, nns)
    counts = zeros(Int, length(rlat))
    _add_counts!(counts, nns)
    itr_Is = findall(==(1), counts)
    all_Is = itr_Is
    iter_count = 1
    while iter_count < maxdepth
        for i in itr_Is
            for j in _iter_neighbor_indices(rlat, nns, i)
                if j != 0 && !(j in all_Is)
                    counts[j] -= 1
                end
            end
        end
        counts[itr_Is] .= 0
        itr_Is = findall(==(1), counts)
        isempty(itr_Is) && break
        mergesorted!(all_Is, itr_Is)
        iter_count += 1
    end
    deleteat!(rlat, all_Is)
end

"""
    removedangling!(lat[, maxdepth])

Remove dangling sites from the lattice. A site is considered dangling if it has less than 2
neighbors. The function will remove all dangling sites and their neighbors recursively up to
`maxdepth` levels — the default is `Inf`.
"""
function removedangling!(lat::AbstractLattice, maxdepth=Inf)
    maxdepth ≤ 0 && return lat
    rlat = addlookuptable(lat)
    _removedangling!(rlat, maxdepth, NearestNeighbor(rlat, 1))
    return lat
end
