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

---
    PointFlux(flux, [point; gauge])

Construct a `PointFlux` object with given flux and point.

The optional `gauge` argument can be used to specify the gauge of the field. Possible values
are `:axial` (``A(r) = B \\times \\frac{r}{|r|}``) and `:singular` (the the phase changes if the
particle passes below the point). The default is `:axial`.
"""
struct PointFlux{GaugeT} <: AbstractField
    flux::Float64
    point::NTuple{2,Float64}
    function PointFlux(Phi,point=(0,0); gauge=:axial)
        gauge in DELTA_GAUGES || throw(ArgumentError("Invalid gauge: $gauge; expected one of $DELTA_GAUGES"))
        new{gauge}(Phi, point)
    end
end
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
