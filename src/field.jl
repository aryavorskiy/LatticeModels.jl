import Base: show, +, *
using LinearAlgebra, StaticArrays, Logging
abstract type AbstractField end

show(io::IO, ::MIME"text/plain", ::T) where {T<:AbstractField} = print(io, "$T field")
@doc raw"""
    vector_potential(field, point)

Returns vector potential $\overrightarrow{A}$ for `field` in location `point`.

This function should be defined for new field types, but it is not necessary
unless you want to use built-in trapezoidal rule integrating.
"""
vector_potential(::FT, point) where {FT<:AbstractField} =
    error("no vector potential function defined for field type $(FT)")

dot_assuming_zeros(m::SVector{M}, n::SVector{N}) where {M,N} = m[1:min(M, N)]' * n[1:min(M, N)]

@doc raw"""
    path_integral(field, p1, p2; n_steps=1)

Calculates the $\int_{p1}^{p2} \overrightarrow{A} \cdot \overrightarrow{dl}$ integral using the trapezoidal rule.
Increase `n_steps` to improve accuracy (for linear fields like Landau or symmetrical calibrations the formula is accurate).
If needed, redefine this function for specific field types - this is likely to boost accuracy and performance.
"""
function path_integral(field::AbstractField, p1, p2; n_steps=1)
    integral = 0.0
    dp = (p2 - p1) / n_steps
    p = p1 + 0.5dp
    for _ in 1:n_steps
        integral += dot_assuming_zeros(dp, SVector(vector_potential(field, p)))
        p += dp
    end
    integral
end

"""
    apply_field!(hamiltonian, field[; nsteps])

Applies magnetic field to given hamiltonian matrix by adjusting the phase factors.
"""
function apply_field!(ham::LatticeOperator, field::AbstractField)
    l = lattice(ham)
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            if i > j && !iszero(ham[i, j])
                p1 = site_coords(l, site1)
                p2 = site_coords(l, site2)
                pmod = exp(2π * im * path_integral(field, p1, p2))
                !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
                ham[i, j] *= pmod
                ham[j, i] *= pmod'
            end
            j += 1
        end
        i += 1
    end
end

function _wrap_block!(f::Function, block::Expr, fields::Vector)
    _begin = 1 + Meta.isexpr(block, :call)
    for i in _begin:length(block.args)
        if block.args[i] isa Symbol && block.args[i] in fields
            block.args[i] = f(block.args[i])
        elseif block.args[i] isa Expr
            _wrap_block!(f, block.args[i], fields)
        end
    end
end

function _extract_varname(var)
    if var isa Symbol
        return var
    elseif Meta.isexpr(var, :(::)) && var.args[1] isa Symbol
        return var.args[1]
    else
        error("cannot extract variable name; invalid definition")
    end
end

"""
    @field_def block

Defines a new magnetic field type.
"""
macro field_def(struct_block)
    if !Meta.isexpr(struct_block, :struct)
        error("Struct block expected")
    end
    struct_head = struct_block.args[2]
    if Meta.isexpr(struct_head, :call)
        struct_name, struct_args... = struct_head.args
    elseif struct_head isa Symbol
        struct_name = struct_head
        struct_args = []
    else
        error("Invalid struct name")
    end
    struct_params = _extract_varname.(struct_args)
    body = struct_block.args[3]
    struct_definition = quote
        import LatticeModels: path_integral, vector_potential
        import Base: show
        struct $struct_name <: AbstractField
            $(struct_args...)
        end
    end
    for statement in body.args
        if Meta.isexpr(statement, (:function, :(=)))
            fn_def, fn_body = statement.args
            _wrap_block!(fn_body, struct_params) do field
                :(field.$field)
            end
            if !Meta.isexpr(fn_def, :call)
                error("not a function definition")
            end
            if length(fn_def.args) < 2
                error("function definition with no arguments")
            end
            fn_name, fn_args... = fn_def.args
            if fn_name === :vector_potential
                if all(isa.(fn_args, Symbol))
                    push!(struct_definition.args, :(
                        function LatticeModels.vector_potential(field::$struct_name, point)
                            $(fn_args...), = point
                            $fn_body
                        end
                    ))
                elseif length(fn_args) == 1 && Meta.isexpr(fn_args[1], :...)
                    fn_arg = fn_args[1].args[1]
                    push!(struct_definition.args, :(
                        function LatticeModels.vector_potential(field::$struct_name, $fn_arg)
                            $fn_body
                        end
                    ))
                else
                    error("invalid argument types")
                end
            elseif fn_name === :path_integral
                local _begin = 1
                Meta.isexpr(fn_args[1], :parameters) && (_begin = 2)
                le = length(fn_args) - _begin + 1
                le != 2 &&
                    error("path_integral must have 2 positional arguments")
                insert!(fn_args, _begin, :(field::$struct_name))
                push!(struct_definition.args, :(
                    function LatticeModels.path_integral($(fn_args...))
                        local vector_potential(p) = LatticeValues.SVector(vector_potential(field, p))
                        $fn_body
                    end
                ))
            elseif fn_name === :show
                push!(struct_definition.args, :(import Base: show),
                    :(function Base.show($(fn_args...), field::$struct_name)
                        $fn_body
                    end
                    ))
            else
                @warn "function definition $fn_def ignored"
            end
        elseif Meta.isexpr(statement, :(:=), 2)
            key, arg = statement.args
            if key === :n_steps
                push!(struct_definition.args, :(
                    LatticeModels.path_integral(field::$struct_name, p1, p2) =
                        LatticeModels.path_integral(field, p1, p2; n_steps=$arg)
                ))
            else
                @warn "unsupported key $key ignored"
            end
        elseif statement isa Expr
            error("not a function definition or key assignment:\n$statement")
        end
    end
    esc(struct_definition)
end

@field_def struct NoField
    path_integral(p1, p2) = 0
end
"""
    NoField <: AbstractField

A stub object representing zero magnetic field.
Use it as a default magnetic field argument in functions - this will not cause any performance overhead.
"""
NoField

@field_def struct LandauField(B::Number)
    vector_potential(x) = SA[0, x*B]
    path_integral(p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * B
    show(io::IO, ::MIME"text/plain") = print(io, "Landau calibration field; B = $B flux quanta per 1×1 plaquette")
end
"""
    LandauField <: AbstractField

An object representing Landau calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
LandauField

@field_def struct SymmetricField(B::Number)
    vector_potential(x, y) = SA[-y, x] * B / 2
    path_integral(p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * B
    show(io::IO, ::MIME"text/plain") = print(io, "Symmetric calibration field; B = $B flux quanta per 1×1 plaquette")
end
"""
    SymmetricField <: AbstractField

An object representing symmetrically calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
SymmetricField

_angle(p1, p2) = asin((1.0 - 1e-11) * det(hcat(p1, p2)) / norm(p1) / norm(p2))
@field_def struct FluxField(B::Number, P::NTuple{2,Number})
    function vector_potential(x, y)
        norm = (x^2 + y^2)
        (-y / norm * B, x / norm * B)
    end
    function path_integral(p1, p2)
        Pv = SVector(P)
        p1 = p1[1:2] - Pv
        p2 = p2[1:2] - Pv
        if iszero(p1) || iszero(p2)
            return 0.0
        end
        _angle(p1, p2) * B
    end
    show(io::IO, ::MIME"text/plain") = print(io, "Delta flux field through point $P; B = $B flux quanta")
end
"""
    FluxField <: AbstractField

An object representing a small magnetic flux through given point. The field is directed along z-axis.
Fields:
- `B`: The magnetic field value
- `point`: A `NTuple{2, Number}` representing the point where the magnetic flux is located.
"""
FluxField

struct FieldSum{N} <: AbstractField
    fields::NTuple{N,AbstractField}
end
vector_potential(f::FieldSum, p1) = sum(SVector(vector_potential(field, p1)) for field in f.fields)
path_integral(f::FieldSum, p1, p2) = sum(path_integral(field, p1, p2) for field in f.fields)
function show(io::IO, m::MIME"text/plain", f::FieldSum{N}) where {N}
    print(io, "Sum of $N fields:\n")
    i = 1
    for field in f.fields
        print("#$i: ")
        show(io, m, field)
        println()
        i += 1
    end
end

+(f::Vararg{AbstractField,N}) where {N} = FieldSum{N}(f)
