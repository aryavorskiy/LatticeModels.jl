import Base: show, +, *
using LinearAlgebra, StaticArrays, Logging
abstract type AbstractField end

show(io::IO, ::MIME"text/plain", ::T) where {T<:AbstractField} = print(io, "$T field")
vector_potential(::FT, point::Vector) where {FT<:AbstractField} =
    error("no vector potential function defined for field type $(MT)")

dot_assuming_zeros(m::SVector{M}, n::SVector{N}) where {M, N} = m[1:min(M, N)]' * n[1:min(M, N)]

function trip_integral(field::AbstractField, p1, p2; n_integrate=1)
    integral = 0.0
    dp = (p2 - p1) / n_integrate
    p = p1 + 0.5dp
    for _ in 1:n_integrate
        integral += dot_assuming_zeros(dp, SVector(vector_potential(field, p)))
        p += dp
    end
    integral
end

function apply_field!(lo::LatticeOperator, field::AbstractField)
    l = lo.basis.lattice
    N = dims_internal(lo.basis)
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            if i > j && !iszero(lo.operator[N*(i-1)+1:N*i, N*(j-1)+1:N*j])
                p1 = coords(l, site1)
                p2 = coords(l, site2)
                pmod = exp(2π * im * trip_integral(field, p1, p2))
                !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
                lo.operator[N*(i-1)+1:N*i, N*(j-1)+1:N*j] *= pmod
                lo.operator[N*(j-1)+1:N*j, N*(i-1)+1:N*i] *= pmod'
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
        import LatticeModels: trip_integral, vector_potential
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
            elseif fn_name === :trip_integral
                local _begin = 1
                Meta.isexpr(fn_args[1], :parameters) && (_begin = 2)
                le = length(fn_args) - _begin + 1
                le != 2 &&
                    error("trip_integral must have 2 positional arguments")
                insert!(fn_args, _begin, :(field::$struct_name))
                push!(struct_definition.args, :(
                    function LatticeModels.trip_integral($(fn_args...))
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
            if key === :n_integrate
                push!(struct_definition.args, :(
                    LatticeModels.trip_integral(field::$struct_name, p1, p2) =
                        LatticeModels.trip_integral(field, p1, p2; n_integrate=$arg)
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
    trip_integral(p1, p2) = 0
end

@field_def struct LandauField(B::Number)
    vector_potential(x) = SA[0, x * B]
    trip_integral(p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * B
    show(io::IO, ::MIME"text/plain") = print(io, "Landau calibration field; B = $B flux quanta per 1×1 plaquette")
end

@field_def struct SymmetricField(B::Number)
    vector_potential(x, y) = SA[-y, x] * B / 2
    trip_integral(p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * B
    show(io::IO, ::MIME"text/plain") = print(io, "Symmetric calibration field; B = $B flux quanta per 1×1 plaquette")
end

_angle(p1, p2) = asin(det(hcat(p1, p2)) / norm(p1) / norm(p2) / (1.0 + 1e-11))
@field_def struct FluxField(B::Number, P::NTuple{2,Number})
    function trip_integral(p1, p2)
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

struct FieldSum{N}
    fields::NTuple{N,AbstractField}
end
trip_integral(f::FieldSum, p1, p2) = sum(trip_integral(field, p1, p2) for field in f.fields)
function show(io::IO, m::MIME"text/plain", f::FieldSum{N}) where {N}
    print(io, "Sum of $N fields:\n")
    i = 1
    for field in f.fields
        print("#$i: ")
        show(io, m, field)
        println()
    end
end

+(f::Vararg{AbstractField,N}) where {N} = FieldSum{N}(f)
