import Base: +, *
using LinearAlgebra, StaticArrays, Logging
abstract type AbstractField end
adapt_field(fld::AbstractField, ::AbstractLattice) = fld

"""
    vector_potential(field, point)

Returns vector potential ``\\overrightarrow{A}`` for `field` in location `point`.

This function should be defined for new field types, but it is not necessary
unless you want to use built-in trapezoidal rule integrating.
"""
vector_potential(fld::AbstractField, _) =
    throw(ArgumentError("no vector potential function defined for field type $(typeof(fld))"))

"""
    line_integral(field, p1, p2[, n_steps=1])

Calculates the ``\\int_{p1}^{p2} \\overrightarrow{A} \\cdot \\overrightarrow{dl}`` integral using the trapezoidal rule.
Increase `n_steps` to improve accuracy (note that for linear field gauges like Landau or symmetrical the formula is already pefrectly accurate).
If needed, redefine this function for specific field types — this is likely to boost accuracy and performance.
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
Use it as a default magnetic field argument in functions — this will not cause any performance overhead.
"""
struct NoField <: AbstractField end
line_integral(::NoField, p1, p2) = 0

"""
    GaugeField <: AbstractField

A gauge field defined by a vector potential function.

---
    GaugeField(func; n)

Create a gauge field with a given vector potential function `func`.

## Arguments
- `func`: a function that takes a point in space and returns the vector potential at this point as a `SVector` or `Tuple`.
- `n`: the number of steps to use in the trapezoidal rule integration.

## Example
```julia
field = GaugeField(n = 10) do p
    (-0.5 * p[2], 0.5 * p[1], 0)
end
```
"""
struct GaugeField{FuncT<:Function} <: AbstractField
    func::FuncT
    n::Int
    GaugeField(func::FuncT; n::Int) where FuncT = new{FuncT}(func, n)
end
vector_potential(field::GaugeField, p1::SVector) = field.func(p1)
line_integral(field::GaugeField, p1, p2) =
    line_integral(field, p1, p2, field.n)

"""
    LineIntegralGaugeField <: AbstractField

A gauge field defined by a line integral function.

---
    LineIntegralGaugeField(func)

Create a gauge field with a given line integral function `func`. The function should take
two points in space and return the line integral of the vector potential between these points.

## Example
```julia
field = LineIntegralGaugeField() do p1, p2
    0.5 * (p1[1] * p2[2] - p1[2] * p2[1])   # A = [-y/2, x/2, 0]
end
```
"""
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
function Base.show(io::IO, f::FieldSum)
    i = 1
    for field in f.fields
        show(io, field)
        i += 1
        i <= length(f.fields) && print(io, " + ")
    end
end

_fldterms(f::FieldSum) = f.fields
_fldterms(f::AbstractField) = (f,)
function Base.add_sum(f1::FieldSum{<:NTuple{32}}, f2::AbstractField)
    @warn "FieldSum is not optimized for a large number of fields; expect compiler issues"
    return f1 + f2
end
+(f1::AbstractField, f2::AbstractField) = FieldSum((_fldterms(f1)..., _fldterms(f2)...))
adapt_field(f::FieldSum, l::AbstractLattice) = FieldSum(adapt_field.(f.fields, Ref(l)))
