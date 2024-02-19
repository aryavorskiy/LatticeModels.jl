import Base: +, *
using LinearAlgebra, StaticArrays, Logging
abstract type AbstractField end

@doc raw"""
    vector_potential(field, point)

Returns vector potential $\overrightarrow{A}$ for `field` in location `point`.

This function should be defined for new field types, but it is not necessary
unless you want to use built-in trapezoidal rule integrating.
"""
vector_potential(fld::AbstractField, _) =
    throw(ArgumentError("no vector potential function defined for field type $(typeof(fld))"))

@doc raw"""
    line_integral(field, p1, p2[, n_steps=1])

Calculates the $\int_{p1}^{p2} \overrightarrow{A} \cdot \overrightarrow{dl}$ integral using the trapezoidal rule.
Increase `n_steps` to improve accuracy (note that for linear field gauges like Landau or symmetrical the formula is already pefrectly accurate).
If needed, redefine this function for specific field types - this is likely to boost accuracy and performance.
"""
line_integral(field::AbstractField, p1, p2) =
    dot_assuming_zeros(p2 - p1, vector_potential(field, (p2 + p1) / 2))
function line_integral(field::AbstractField, p1, p2, n_steps)
    integral = 0.0
    dp = (p2 - p1) / n_steps
    p = p1 + 0.5dp
    for _ in 1:n_steps
        integral += dot_assuming_zeros(dp, vector_potential(field, p))
        p += dp
    end
    integral
end

"""
    NoField <: AbstractField

A stub object representing zero magnetic field.
Use it as a default magnetic field argument in functions - this will not cause any performance overhead.
"""
struct NoField <: AbstractField end
line_integral(::NoField, p1, p2) = 0

struct GaugeField{FuncT<:Function} <: AbstractField
    func::FuncT
    n::Int
    GaugeField(func::FuncT; n::Int) where FuncT = new{FuncT}(func, n)
end
vector_potential(field::GaugeField, p1::SVector) = field.func(p1)
line_integral(field::GaugeField, p1, p2) =
    line_integral(field, p1, p2, field.n)

struct LineIntegralGaugeField{FuncT<:Function} <: AbstractField
    func::FuncT
    LineIntegralGaugeField(func::FuncT) where FuncT = new{FuncT}(func)
end
line_integral(field::LineIntegralGaugeField, p1, p2) = field.func(p1, p2)

struct FieldSum{FT<:Tuple} <: AbstractField
    fields::FT
end
vector_potential(f::FieldSum, p1) = sum(SVector(vector_potential(field, p1)) for field in f.fields)
line_integral(f::FieldSum, p1, p2) = sum(line_integral(field, p1, p2) for field in f.fields)
function Base.show(io::IO, m::MIME"text/plain", f::FieldSum)
    print(io, "Sum of ", format_number(length(f.fields), "gauge field term"))
    i = 1
    for field in f.fields
        print(io, "\n#$i: ")
        show(io, m, field)
        i += 1
    end
end

+(f::Vararg{AbstractField}) = FieldSum{typeof(f)}(f)
