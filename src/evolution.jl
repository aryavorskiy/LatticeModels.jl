import Base: exp
using LinearAlgebra

function taylor_exp(A::AbstractMatrix, k::Int)
    B = one(A) + A
    if k == 1
        return B
    end
    M = copy(A)
    for i in 2:k
        M *= A / i
        B += M
    end
    return B
end

exp(A::LatticeOperator) = LatticeArray(A.basis, exp(A.operator))
taylor_exp(A::LatticeOperator, k::Int) = LatticeArray(A.basis, taylor_exp(A.operator, k))

@doc raw"""
    evolution_operator(H, t[, k])

Calculates the unitary evolution operator using the formula

$ \mathcal{U}(t) = e^{-\frac{1}{i\hbar} \hat{H} t} $

# Arguments
- `H`: the hamiltonian matrix
- `t`: the evolution time
- `k`: if provided, the exponent will be calculated using a Taylor series expansion with order k
"""
evolution_operator(H, t::Real) = exp((-im * t) * H)
evolution_operator(H, t::Real, k::Int) = taylor_exp((-im * t) * H, k)

evolved(P::AbstractMatrix, ev::AbstractMatrix) = ev * P * ev'
evolved(V::AbstractVector, ev::AbstractMatrix) = ev * V
evolved(LA, ev) = @on_lattice evolved(LA, ev)

function _expand_chain(chain)
    if Meta.isexpr(chain, :-->)
        (_expand_chain(chain.args[1])..., _expand_chain(chain.args[2])...)
    else
        (chain,)
    end
end

_expr_depends_on(::Any, ::Symbol) = true
_expr_depends_on(expr::Symbol, dep_var::Symbol) = (expr === dep_var)
function _expr_depends_on(expr::Expr, dep_var::Symbol)
    _begin = 1 + Meta.isexpr(expr, (:call, :->))
    any(_expr_depends_on(e, dep_var) for e in expr.args[_begin:end])
end

function _evolution_operator_call(H_sym, dt_sym, k)
    if k === nothing
        :(evolution_operator($H_sym, $dt_sym))
    else
        :(evolution_operator($H_sym, $dt_sym, $k))
    end
end

function _evolution_block(rules, loop; k=nothing, rtol=1e-12)
    if !Meta.isexpr(loop, :for)
        error("for loop expected")
    end
    loop_iter, loop_body = loop.args
    if !(Meta.isexpr(loop_iter, :(=)))
        error("invalid iteration specification '$loop_iter'")
    end
    loop_var, loop_range = loop_iter.args
    !Meta.isexpr(loop_body, :block) && error("malformed loop body")

    inits = []
    ham_evals = []
    evolutor_updates = []
    p_evolutions = []

    if typeof(rules) != Expr || rules.head ∉ (:braces, :bracescat)
        error("evolution specifier list should be a braces notation, not '$(loop.head)'")
    end

    hamiltonian_functions = []
    hamiltonian_aliases = Dict{Symbol,Int}()
    for statement in rules.args
        if Meta.isexpr(statement, :(:=), 2)
            ham_sym, ham_expr = statement.args
            !(ham_sym isa Symbol) && error("assignment lvalue must be a Symbol")
            if ham_sym in keys(hamiltonian_aliases)
                error("redefinition of alias $statement not allowed
                (previous $ham_sym := $(hamiltonian_functions[hamiltonian_aliases[ham_sym]]))")
            end
            !(ham_expr in hamiltonian_functions) && push!(hamiltonian_functions, ham_expr)
            hamiltonian_aliases[ham_sym] = findfirst(==(ham_expr), hamiltonian_functions)
        elseif statement isa Expr
            chain = _expand_chain(statement)
            if length(chain) == 3
                !(chain[2] in keys(hamiltonian_aliases)) &&
                    push!(hamiltonian_functions, chain[2])
            else
                error("invalid specifier:\n$statement")
            end
        end
    end
    ham_i = 1
    for ham_expr in hamiltonian_functions
        h_eval = Symbol("ham_eval_$ham_i")
        h_eval_new = Symbol("ham_eval_new_$ham_i")
        p_target_ev = Symbol("evolutor_$ham_i")
        if _expr_depends_on(ham_expr, loop_var)
            push!(inits,
                :(local $h_eval = $(esc(ham_expr))),
                :(local $h_eval_new = $(esc(ham_expr))),
                :(local $p_target_ev = _unwrap_from_macro(one, $(esc(ham_expr)))))
            push!(ham_evals, :($h_eval_new = $(esc(ham_expr))))
            push!(evolutor_updates,
                :(
                    if dt_changed || $h_eval != $h_eval_new
                        $p_target_ev = $(_evolution_operator_call(h_eval_new, :dt, k))
                    end
                ),
                :($h_eval = $h_eval_new))
        else
            push!(inits,
                :(local $h_eval = $(esc(ham_expr))),
                :(local $p_target_ev = _unwrap_from_macro(one, $(esc(ham_expr)))))
            push!(evolutor_updates,
                :(
                    if dt_changed
                        $p_target_ev = $(_evolution_operator_call(h_eval, :dt, k))
                    end
                ))
        end
        ham_i += 1
    end

    for statement in rules.args
        if statement isa LineNumberNode
            continue
        elseif Meta.isexpr(statement, :(:=), 2)
            local ham_var, ham_expr = statement.args
            ham_i = findfirst(==(ham_expr), hamiltonian_functions)
            h_eval = Symbol("ham_eval_$ham_i")
            h_eval_new = Symbol("ham_eval_new_$ham_i")
            if _expr_depends_on(ham_expr, loop_var)
                push!(ham_evals, :($(esc(ham_var)) = $h_eval_new))
            else
                push!(ham_evals, :($(esc(ham_var)) = $h_eval))
            end
        elseif statement isa Expr
            p_initial, ham_expr, p_target = _expand_chain(statement)
            ham_i = get(hamiltonian_aliases, ham_expr, findfirst(==(ham_expr), hamiltonian_functions))
            ham_i === nothing && error()
            h_eval = Symbol("ham_eval_$ham_i")
            h_eval_new = Symbol("ham_eval_new_$ham_i")
            p_target_ev = Symbol("evolutor_$ham_i")
            push!(inits,
                :(local $(esc(p_target)) = _unwrap_from_macro(copy, $(esc(p_initial)))))
            push!(p_evolutions,
                :($(esc(p_target)) =
                    evolved($(esc(p_target)), $p_target_ev)))
        end
    end
    quote
        local $(esc(loop_var)) = first($(esc(loop_range)))
        $(inits...)
        local t_inner = zero(eltype($(esc(loop_range))))
        local dt_old = zero(eltype($(esc(loop_range))))
        local dt_changed::Bool = false
        for $(esc(loop_var)) in $(esc(loop_range))
            local dt = $(esc(loop_var)) - t_inner
            dt == 0 && continue
            dt_changed = abs(dt - dt_old) / dt > $rtol
            if dt_changed
                dt_old = dt
            end
            $(ham_evals...)
            $(evolutor_updates...)
            $(p_evolutions...)
            $(esc(loop_body))
            t_inner = $(esc(loop_var))
        end
    end
end

"""
    @evolution [kwargs...] {rules...} for_loop

Generates an environment with defined hamiltonian and density matrices that evolve by certain laws.
See [Unitary evolution](evolution.md) for more details.
"""
macro evolution(args...)
    rules, loop = args[end-1:end]
    kwargs = Dict(a.args for a in args[begin:end-2])
    _evolution_block(rules, loop; kwargs...)
end
