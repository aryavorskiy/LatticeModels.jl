import Base: +

"""
    LandauGauge <: AbstractField

An object representing Landau gauge uniform magnetic field along z-axis.

## Fields
- `B`: The magnetic field value
"""
struct LandauGauge <: AbstractField
    B::Float64
end
vector_potential(field::LandauGauge, p1) = (0, p1[1] * field.B)
line_integral(field::LandauGauge, p1, p2) = (p1[1] + p2[1]) * (p2[2] - p1[2]) * field.B / 2
Base.show(io::IO, ::MIME"text/plain", field::LandauGauge) =
    print(io, "Landau gauge uniform field; B = $(field.B) flux quanta per 1×1 area")

"""
    SymmetricGauge <: AbstractField

An object representing symmetrically gauged uniform magnetic field along z-axis.

## Fields
- `B`: The magnetic field value
"""
struct SymmetricGauge <: AbstractField
    B::Float64
end
vector_potential(field::SymmetricGauge, p1) = SA[-p1[2], p1[1]] * field.B / 2
line_integral(field::SymmetricGauge, p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * field.B
Base.show(io::IO, ::MIME"text/plain", field::SymmetricGauge) =
    print(io, "Symmetric gauge uniform field; B = $(field.B) flux quanta per 1×1 area")

const DELTA_GAUGES = (:axial, :singular)
"""
    PointFlux{GaugeT} <: AbstractField

An object representing a small magnetic flux through given point. The field is directed along z-axis.


## Fields
- `flux`: The magnetic flux value.
- `point`: A `Tuple` of x and y coordinates of the point.
"""
struct PointFlux{GaugeT} <: AbstractField
    flux::Float64
    point::NTuple{2,Float64}
    function PointFlux{GaugeT}(phi, point=(0,0)) where GaugeT
        GaugeT in DELTA_GAUGES || throw(ArgumentError("Invalid gauge: $GaugeT; expected one of $DELTA_GAUGES"))
        new{GaugeT}(phi, point)
    end
end

"""
    PointFlux(flux, [point; gauge])

Construct a `PointFlux` object with given flux and point.

The optional `gauge` argument can be used to specify the gauge of the field. Possible values
are `:axial` (``A(r) = B \\times \\frac{r}{|r|}``) and `:singular` (the the phase changes if the
particle passes below the point). The default is `:axial`.
"""
PointFlux(phi, point=(0,0); gauge=:axial) = PointFlux{gauge}(phi, point)

function vector_potential(field::PointFlux{:axial}, p1)
    (x, y) = p1
    normsq = (x^2 + y^2)
    SA[-y, x] / normsq / 2pi * field.flux
end
const FLUX_WARNING = "Hopping line goes through the point of the delta flux, numerical errors possible"
function line_integral(field::PointFlux{:axial}, p1, p2)
    Pv = SVector(field.point)
    p1 = p1[1:2] - Pv
    p2 = p2[1:2] - Pv
    if norm(p1) < 1e-11 || norm(p2) < 1e-11
        @warn FLUX_WARNING
        return zero(field.flux)
    end
    nnorm = norm(p1) * norm(p2)
    anglesinsign = det(hcat(p1, p2))
    anglecos = dot(p1, p2) / nnorm / (1 + 1e-11)
    angle = acos(anglecos) * sign(anglesinsign)
    return angle * field.flux / 2pi
end
function vector_potential(::PointFlux{:singular}, p1)
    throw(ArgumentError("Vector potential for a point flux is singular"))
end
function line_integral(field::PointFlux{:singular}, p1, p2)
    Pv = SVector(field.point)
    x1, y1 = p1[1:2] - Pv
    x2, y2 = p2[1:2] - Pv
    sg = x2 - x1
    if abs(sg) < 1e-11
        abs(x1) < 1e-11 && y1 * y2 ≤ 0 && @warn FLUX_WARNING
        return zero(field.flux)
    end
    if x1 * x2 > 0 || iszero(max(x1, x2))
        return zero(field.flux)
    end
    yintercept = (- y1 * x2 + y2 * x1) / (x1 - x2)
    abs(yintercept) < 1e-11 && @warn FLUX_WARNING
    yintercept > 0 ? zero(field.flux) : field.flux * sign(sg)
end
Base.show(io::IO, field::PointFlux{GaugeT}) where GaugeT =
    print(io, "PointFlux(", field.flux, ", ", field.point, "; gauge=:", GaugeT, ")")
Base.show(io::IO, ::MIME"text/plain", field::PointFlux{GaugeT}) where GaugeT =
    print(io, "Point flux field through point $(field.point), $(GaugeT) gauge; Φ = $(field.flux) flux quanta")

"""
    PointFluxes{GaugeT} <: AbstractField

An object representing a collection of small magnetic fluxes through given points. The field is directed along z-axis.

## Fields
- `fluxes`: The magnetic flux values.
- `points`: A vector of `NTuple{2,Float64}` points.
"""
struct PointFluxes{GaugeT} <: AbstractField
    fluxes::Vector{Float64}
    points::Vector{NTuple{2,Float64}}
    function PointFluxes{GaugeT}(fluxes::AbstractVector, points::AbstractVector) where GaugeT
        GaugeT in DELTA_GAUGES || throw(ArgumentError("Invalid gauge: $GaugeT; expected one of $DELTA_GAUGES"))
        length(fluxes) == length(points) || throw(ArgumentError("Length of fluxes and points should be the same"))
        new{GaugeT}(fluxes, points)
    end
end

"""
    PointFluxes([fluxes, points; offset=(0, 0), gauge=:axial])

Construct a `PointFluxes` object with given fluxes and points.

The optional `gauge` argument can be used to specify the gauge of the field.

## Arguments
- `fluxes`: A vector of flux values. Also can be a single value, in which case it will be broadcasted to all points.
- `points`: A vector of points or a `AbstractLattice` object. In the latter case the sites will be interpreted as points.

If both arguments are omitted, an empty field is created. You can add fluxes to it later using `push!` or `append!`.

## Keyword arguments
- `offset`: An offset to add to the points, default is `(0, 0)`. Valid only if `points` is a `AbstractLattice`, otherwise an error is thrown.
- `gauge`: The gauge of the field. Possible values are `:axial` and `:singular`. The default is `:axial`.
"""
PointFluxes(fluxes::AbstractVector, points::AbstractVector; gauge=:axial) =
    PointFluxes{gauge}(fluxes, points)
PointFluxes(flux::Real, points::AbstractVector; kw...) =
    PointFluxes(fill(Float64(flux), length(points)), points; kw...)
PointFluxes(fluxes::Any, points::AbstractLattice; offset=(0, 0), kw...) =
    return PointFluxes(fluxes, [(x, y) .+ offset for (x, y) in points]; kw...)
PointFluxes(;gauge=:axial) = PointFluxes{gauge}([], [])

_fluxgauge(::PointFlux{GaugeT}) where GaugeT = GaugeT
_fluxgauge(::PointFluxes{GaugeT}) where GaugeT = GaugeT
_fluxgauge(field::Any) = throw(TypeError(:PointFluxes, "_fluxgauge", PointFlux, field))
"""
    PointFluxes(fields[; gauge])

Construct a `PointFluxes` object from a collection of `PointFlux` objects.

An error is thrown if the gauges of the fields are inconsistent. You can specify the gauge
of the field explicitly using the `gauge` keyword argument.

See also: [`PointFlux`](@ref).
"""
function PointFluxes(fields; gauge=nothing)
    isempty(fields) && throw(ArgumentError("Empty field collection"))
    force_gauge = gauge === nothing ? _fluxgauge(first(fields)) : gauge
    fluxes = Float64[]
    points = NTuple{2,Float64}[]
    for field in fields
        if gauge === nothing && _fluxgauge(field) !== force_gauge
            throw(ArgumentError("Inconsistent gauges: $gauge vs $force_gauge"))
        end
        push!(fluxes, field.flux)
        push!(points, field.point)
    end
    PointFluxes{force_gauge}(fluxes, points)
end
vector_potential(field::PointFluxes{GaugeT}, p1) where {GaugeT} =
    sum(vector_potential(PointFlux{GaugeT}(flux, point), p1) for (flux, point) in zip(field.fluxes, field.points))
line_integral(field::PointFluxes{GaugeT}, p1, p2) where {GaugeT} =
    sum(line_integral(PointFlux{GaugeT}(flux, point), p1, p2) for (flux, point) in zip(field.fluxes, field.points))
Base.show(io::IO, field::PointFluxes{GaugeT}) where GaugeT =
    print(io, "PointFluxes(", !isempty(field.fluxes) && allequal(field.fluxes) ?
    first(field.fluxes) : field.fluxes, ", ", field.points, "; gauge=:", GaugeT, ")")

"""
    push!(fields, field)

Add a `PointFlux` object to a `PointFluxes` object. An error is thrown if the gauges of the fields are inconsistent.
"""
function Base.push!(fields::PointFluxes, flux::PointFlux, fluxes...)
    if _fluxgauge(flux) !== _fluxgauge(fields)
        throw(ArgumentError("Inconsistent gauges: $(_fluxgauge(flux)) vs $(_fluxgauge(fields))"))
    end
    push!(fields.fluxes, flux.flux)
    push!(fields.points, flux.point)
    isempty(fluxes) || push!(fields, fluxes...)
    return fields
end

"""
    append!(fields, fields2)

Add a `PointFluxes` object to another `PointFluxes` object. An error is thrown if the gauges of the fields are inconsistent.
"""
function Base.append!(fields::PointFluxes, fields2::PointFluxes)
    if _fluxgauge(fields2) !== _fluxgauge(fields)
        throw(ArgumentError("Inconsistent gauges: $(_fluxgauge(fields2)) vs $(_fluxgauge(fields))"))
    end
    append!(fields.fluxes, fields2.fluxes)
    append!(fields.points, fields2.points)
end

"""
    periodic_fluxes(l, fl)

Construct a `PointFluxes` object by periodic replication of a point flux over a Bravais lattice.

## Arguments
- `l`: A Bravais lattice.
- `fl`: A `PointFlux` object.
"""
function periodic_fluxes(l::MaybeWithMetadata{BravaisLattice}, fl::PointFlux)
    uv = unitvectors(l)[:, 1:2]
    lv0 = round.(uv \ SVector(fl.point))
    p0 = SVector(fl.point) - uv * lv0
    points = Set{NTuple{2,Float64}}()
    for site in l
        push!(points, Tuple(p0 + uv * site.latcoords))
    end
    return PointFluxes{_fluxgauge(fl)}(fill(fl.flux, length(points)), collect(points))
end
