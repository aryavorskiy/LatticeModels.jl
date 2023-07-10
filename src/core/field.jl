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
    error("no vector potential function defined for field type $(typeof(fld))")

@doc raw"""
    line_integral(field, p1, p2[, n_steps=1])

Calculates the $\int_{p1}^{p2} \overrightarrow{A} \cdot \overrightarrow{dl}$ integral using the trapezoidal rule.
Increase `n_steps` to improve accuracy (note that for linear fields like Landau or symmetrical calibrations the formula is already pefrectly accurate).
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
        struct_name = esc(struct_name)
    elseif struct_head isa Symbol
        struct_name = esc(struct_head)
        struct_args = []
    else
        error("Invalid struct name")
    end
    struct_fields = map(struct_args) do expr
        Meta.isexpr(expr, (:(=), :kw), 2) ? expr.args[1] : expr
    end
    struct_params = map(struct_fields) do var
        var isa Symbol && return var
        Meta.isexpr(var, :(::)) && var.args[1] isa Symbol && return var.args[1]
        error("cannot extract variable name; invalid definition '$var'")
    end
    body = struct_block.args[3]
    struct_definition = quote
        import LatticeModels: line_integral, vector_potential
        import Base: show
        struct $struct_name <: AbstractField
            $(esc.(struct_fields)...)
            $struct_name($(esc.(struct_args)...)) = new($(esc.(struct_fields)...))
        end
    end
    for statement in body.args
        if Meta.isexpr(statement, (:function, :(=)))
            fn_def, fn_body = statement.args
            _wrap_block!(fn_body, struct_params) do field
                :(field.$field)
            end
            fn_body = esc(fn_body)
            if !Meta.isexpr(fn_def, :call)
                error("not a function definition")
            end
            if length(fn_def.args) < 2
                error("function definition with no arguments")
            end
            fn_name, fn_args... = fn_def.args
            escf = :($(esc(:field))::$struct_name)
            if fn_name === :vector_potential
                if all(isa.(fn_args, Symbol))
                    push!(struct_definition.args, :(
                        function LatticeModels.vector_potential($escf, point)
                            $(esc.(fn_args)...), = point
                            $fn_body
                        end
                    ))
                elseif length(fn_args) == 1 && Meta.isexpr(fn_args[1], :...)
                    fn_arg = esc(fn_args[1].args[1])
                    push!(struct_definition.args, :(
                        function LatticeModels.vector_potential($escf, $fn_arg)
                            $fn_body
                        end
                    ))
                else
                    error("invalid argument types")
                end
            elseif fn_name === :line_integral
                Meta.isexpr(fn_args[1], :parameters) &&
                    error("line_integral must not accept kwargs")
                length(fn_args) != 2 &&
                    error("line_integral must have 2 positional arguments")
                push!(struct_definition.args, :(
                    function LatticeModels.line_integral($escf, $(esc.(fn_args)...))
                        $fn_body
                    end
                ))
            elseif fn_name === :show
                push!(struct_definition.args, :(import Base: show),
                    :(function Base.show($(esc.(fn_args)...), $escf)
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
                    LatticeModels.line_integral(field::$struct_name, p1, p2) =
                        LatticeModels.line_integral(field, p1, p2, $arg)
                ))
            else
                @warn "unsupported key $key ignored"
            end
        elseif statement isa Expr
            error("not a function definition or key assignment:\n$statement")
        end
    end
    struct_definition
end

@field_def struct NoField
    line_integral(p1, p2) = 0
end
"""
    NoField <: AbstractField

A stub object representing zero magnetic field.
Use it as a default magnetic field argument in functions - this will not cause any performance overhead.
"""
NoField

struct FieldSum{FT<:Tuple} <: AbstractField
    fields::FT
end
vector_potential(f::FieldSum, p1) = sum(SVector(vector_potential(field, p1)) for field in f.fields)
line_integral(f::FieldSum, p1, p2) = sum(line_integral(field, p1, p2) for field in f.fields)
function Base.show(io::IO, m::MIME"text/plain", f::FieldSum)
    print(io, "Sum of $(length(f.fields)) fields:")
    i = 1
    for field in f.fields
        print("\n#$i: ")
        show(io, m, field)
        i += 1
    end
end

+(f::Vararg{AbstractField}) = FieldSum{typeof(f)}(f)
