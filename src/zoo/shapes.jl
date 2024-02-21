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

function fill_shapes(uc::UnitCell{Sym,N} where Sym, shapes::AbstractShape...; sites::Nullable{Int}=nothing, kw...) where N
    bps = BravaisPointer{N}[]
    if sites === nothing
        new_shapes = shapes
    else
        n_estimate = sum(volume, shapes) * length(uc) / det(unitvectors(uc))
        scale_factor = (sites / n_estimate) ^ (1 / dims(uc))
        new_shapes = scale.(shapes, scale_factor)
    end
    for shape in new_shapes
        add_bravaispointers!(site -> inshape(shape, site), bps, uc, bounding_region(uc, shape))
    end
    b = BravaisLattice(uc, bps)
    return finalize_lattice(b; kw...)
end
fill_shapes(::Type{LT}, shapes::AbstractShape...; kw...) where LT<:BravaisLattice =
    fill_shapes(construct_unitcell(LT), shapes...; kw...)

(::Type{LT})(shapes::AbstractShape...; kw...) where LT<:BravaisLattice =
    fill_shapes(construct_unitcell(LT), shapes...; kw...)

struct Circle{N} <: AbstractShape{N}
    radius::Float64
    center::SVector{N, Float64}
    Circle(radius::Real, center::AbstractVector{<:Real}=zero(SVector{2})) =
        new{length(center)}(radius, center)
end
Circle(center::AbstractVector{<:Real}=zero(SVector{2})) = Circle(1, center)
function circumscribed_sphere(uc::UnitCell, f::Circle)
    return f.radius, add_assuming_zeros(zero(SVector{dims(uc)}), f.center)
end
scale(f::Circle, c) = Circle(c * f.radius, c * f.center)
volume(f::Circle{N}) where N = ndvol(dims(f)) * f.radius^N
inshape(f::Circle, site) = sum(abs2, add_assuming_zeros(site.coords, -f.center)) ≤ f.radius^2

struct Polygon{N} <: AbstractShape{2}
    radius::Float64
    center::SVector{2, Float64}
    function Polygon{N}(radius::Real, center::AbstractVector{<:Real}=zero(SVector{2})) where N
        @check_size center 2
        new{N}(radius, center)
    end
end
Polygon{N}(center::AbstractVector{<:Real}=zero(SVector{2})) where N = Polygon{N}(1, center)
function circumscribed_sphere(::UnitCell, f::Polygon)
    return f.radius, f.center
end
scale(f::Polygon{N}, c) where N = Polygon{N}(c * f.radius, c * f.center)
volume(f::Polygon{N}) where N = N / 2 * sin(2pi / N) * f.radius^2
@generated function inshape(f::Polygon{N}, site) where N
    exprs = [:(y * $(cos(2pi * i / N)) + x * $(sin(2pi * i / N)) < h) for i in 1:N]
    and_expr = reduce((x, y) -> :($x && $y), exprs)
    quote
        h = f.radius * $(cos(pi / N))
        x, y = site.coords - f.center
        return $and_expr
    end
end

const Triangle = Polygon{3}
const Hexagon = Polygon{6}

"""
    GrapheneRibbon(len, wid[, center; kw...])

Construct a graphene ribbon sample with zigzag edges.
To get armchair edges, simply rotate the lattice by 90 degrees.

## Arguments:
- `len`: the length of the ribbon.
- `wid`: the width of the ribbon.
- `center`: the unit cell coordinates of the bottom-left corner of the ribbon. Default is `(0, 0)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function GrapheneRibbon(len, wid, center=(0, 0); kw...)
    j1ind = center[1] - wid:len + center[1] - 1
    j2ind = center[2]:wid + center[2] - 1
    HoneycombLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2 / 2
        if -0.5 ≤ j_prj ≤ len - 0.5
            return !(site.basindex == 2 && j_prj == len - 0.5 ||
                site.basindex == 1 && j_prj == -0.5)
        else
            return false
        end
    end
end
