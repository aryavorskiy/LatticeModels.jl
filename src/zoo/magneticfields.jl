"""
    LandauField <: AbstractField

An object representing Landau gauge uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
struct LandauField <: AbstractField
    B::Float64
end
vector_potential(field::LandauField, p1) = (0, p1[1] * field.B)
line_integral(field::LandauField, p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * field.B
Base.show(io::IO, ::MIME"text/plain", field::LandauField) = print(io, "Landau gauge uniform field; B = $(field.B) flux quanta per 1×1 plaquette")

"""
    SymmetricField <: AbstractField

An object representing symmetrically gauged uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
struct SymmetricField <: AbstractField
    B::Float64
end
vector_potential(field::SymmetricField, p1) = SA[-p1[2], p1[1]] * field.B / 2
line_integral(field::SymmetricField, p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * field.B
Base.show(io::IO, ::MIME"text/plain", field::SymmetricField) = print(io, "Symmetric gauge uniform field; B = $(field.B) flux quanta per 1×1 plaquette")

"""
    FluxField <: AbstractField

An object representing a small magnetic flux through given point. The field is directed along z-axis.
Fields:
- `B`: The magnetic field value
- `point`: A `NTuple{2, Number}` representing the point where the magnetic flux is located.
"""
struct FluxField <: AbstractField
    B::Float64
    P::NTuple{2,Float64}
    FluxField(B,P=(0,0)) = new(B, P)
end
function vector_potential(field::FluxField, p1)
    (x, y) = p1
    normsq = (x^2 + y^2)
    SA[-y, x] / normsq * field.B
end
function line_integral(field::FluxField, p1, p2)
    Pv = SVector(field.P)
    p1 = p1[1:2] - Pv
    p2 = p2[1:2] - Pv
    if iszero(p1) || iszero(p2)
        return 0.0
    end
    asin((1 - 1e-11) * det(hcat(p1, p2)) / norm(p1) / norm(p2)) * field.B
end
Base.show(io::IO, ::MIME"text/plain", field::FluxField) = print(io, "Delta flux field through point $(field.P); B = $(field.B) flux quanta")
