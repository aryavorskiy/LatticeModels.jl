import Base: exp
using LinearAlgebra, ProgressMeter

function myexp!(P::AbstractMatrix, A::AbstractMatrix, factor; threshold=1e-6, nonzero_tol=1e-14)
    mat_norm = norm(A, Inf)
    (mat_norm ≤ 0 || iszero(factor)) && return one(A)
    scaling_factor = nextpow(2, mat_norm * abs(factor))
    A /= scaling_factor
    delta = one(mat_norm)
    LinearAlgebra.checksquare(A)

    P .= one(A)
    next_term = copy(P)
    nt_buffer = similar(next_term)
    n = 1
    while delta > threshold
        if issparse(A)
            next_term = A * next_term
            next_term .*= factor / n
            droptol!(next_term, nonzero_tol)
        else
            mul!(nt_buffer, A, next_term, factor / n, 0)
            copyto!(next_term, nt_buffer)
        end

        delta = norm(next_term, Inf)
        P .+= next_term
        n += 1
    end
    for _ in 1:log2(scaling_factor)
        if issparse(A)
            P = P * P
        else
            mul!(nt_buffer, P, P)
            copyto!(P, nt_buffer)
        end
    end
    return P
end

@doc raw"""
    evolution_operator(H, t)

Calculates the unitary evolution operator using the formula

$ \mathcal{U}(t) = e^{-\frac{1}{i\hbar} \hat{H} t} $

# Arguments
- `H`: the hamiltonian matrix
- `t`: the evolution time
"""
evolution_operator!(Ev::AbstractMatrix, H::AbstractMatrix, t::Real) = myexp!(Ev, H, -im * t)
evolution_operator!(Ev::DataOperator, H::DataOperator, t) = Operator(basis(H), evolution_operator!(Ev.data, H.data, t))

evolved(ket::Ket, ev) = ev * ket
evolved(op::DataOperator, ev) = ev * op * ev'
evolved(bra::Bra, ev) = evolved(bra', ev)'

function _expand_chain(chain)
    if Meta.isexpr(chain, :-->)
        (_expand_chain(chain.args[1])..., _expand_chain(chain.args[2])...)
    else
        (chain,)
    end
end

_expr_depends_on(::Any, ::Symbol) = true
_expr_depends_on(::Number, ::Symbol) = false
_expr_depends_on(::String, ::Symbol) = false
_expr_depends_on(expr::Symbol, dep_var::Symbol) = (expr === dep_var)
function _expr_depends_on(expr::Expr, dep_var::Symbol)
    _begin = 1 + Meta.isexpr(expr, (:call, :->))
    any(_expr_depends_on(e, dep_var) for e in expr.args[_begin:end])
end

function parse_timesequences!(expr, t, seqnames=[])
    !Meta.isexpr(expr, (:if, :elseif, :for, :while, :block)) && return seqnames
    for i in eachindex(expr.args)
        if Meta.isexpr(expr.args[i], :call) && expr.args[i].args[1] == :<--
            seq, val = expr.args[i].args[2:end]
            ref = findfirst(x -> x.first == seq, seqnames)
            if ref === nothing
                seqbox = gensym("$(seq)_box")
                push!(seqnames, seq => seqbox)
            else seqbox = seqnames[ref].second
            end
            expr.args[i] = :($seqbox[$t] = $val)
        else
            parse_timesequences!(expr.args[i], t, seqnames)
        end
    end
    return seqnames
end

one_keeptype(op::DataOperator) = Operator(basis(op), one(op.data))
one_keeptype(mat) = one(mat)

function _evolution_block(rules, loop; k=nothing, rtol=1e-12, show_progress=true, pade=false)
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

    !Meta.isexpr(rules, (:braces, :bracescat)) &&
        error("evolution specifier list should be a braces notation, not '$(rules.head)'")

    hamiltonian_functions = []
    hamiltonian_aliases = Dict{Symbol,Int}()
    for statement in rules.args
        if Meta.isexpr(statement, :(:=), 2)
            ham_sym, ham_expr = statement.args
            ham_sym isa Symbol || error("hamiltonian alias must be a Symbol")
            if ham_sym in keys(hamiltonian_aliases)
                error("""cannot overwrite alias '$ham_sym' with value '$ham_expr'
                ('$(hamiltonian_functions[hamiltonian_aliases[ham_sym]])' assigned before)""")
            end
            ham_expr in hamiltonian_functions || push!(hamiltonian_functions, ham_expr)
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
                :(local $p_target_ev = (1. + 0im) * one_keeptype($(esc(ham_expr)))))
            push!(ham_evals, :($h_eval_new = $(esc(ham_expr))))
            push!(evolutor_updates,
                :(
                    if dt_changed || $h_eval != $h_eval_new
                        $p_target_ev = evolution_operator!($p_target_ev, $h_eval_new, dt)
                    end
                ),
                :($h_eval = $h_eval_new))
        else
            push!(inits,
                :(local $h_eval = $(esc(ham_expr))),
                :(local $p_target_ev = (1. + 0im) * one_keeptype($(esc(ham_expr)))))
            push!(evolutor_updates,
                :(
                    if dt_changed
                        $p_target_ev = evolution_operator!($p_target_ev, $h_eval, dt)
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
                :(local $(esc(p_target)) = copy($(esc(p_initial)))))
            push!(p_evolutions,
                :($(esc(p_target)) =
                    evolved($(esc(p_target)), $p_target_ev)))
        end
    end
    seqs = parse_timesequences!(loop_body, loop_var)
    afters = []
    for (seq, seqbox) in seqs
        push!(inits,
        :(local $(esc(seqbox)) = TimeSequenceContainer()))
        push!(afters,
        :($(esc(seq)) = $(esc(seqbox)).seq))
    end
    quote
        ProgressMeter.ijulia_behavior(:clear)
        local $(esc(loop_var)) = first($(esc(loop_range)))
        $(inits...)
        local t_inner = zero(eltype($(esc(loop_range))))
        local dt_old = zero(eltype($(esc(loop_range))))
        local dt_changed::Bool = false
        local p = Progress(length($(esc(loop_range))), desc="Evolution... ",
            barglyphs=BarGlyphs("[=> ]"), showspeed=true, enabled=$show_progress)
        local tstart = time()
        local dt_evol = 0.
        for $(esc(loop_var)) in $(esc(loop_range))
            tstartevol = time()
            local dt = $(esc(loop_var)) - t_inner
            dt_changed = abs(dt - dt_old) / dt > $rtol
            if dt_changed
                dt_old = dt
            end
            $(ham_evals...)
            $(evolutor_updates...)
            $(p_evolutions...)
            dt_evol += time() - tstartevol
            $(esc(loop_body))
            ProgressMeter.next!(p; showvalues = [("% of time performing evolution",
                round(100 * dt_evol / (time() - tstart), digits=1))])
            t_inner = $(esc(loop_var))
        end
        $(afters...)
    end
end

"""
    @evolution [kwargs...] {rules...} for_loop

Generates an environment with defined hamiltonian and density matrices that evolve by certain laws.
See [Unitary evolution](evolution.md) for more details.

**Keyword arguments:**
- `k`: order of the Taylor expansion for matrix exponent. If omitted, the default `exp` function will be used.
- `pade`: set this to true to use Padé approximant formula instead of Taylor expansion.
- `rtol`: the relative tolerance to decide whether the `Δt` changed between iterations or not. `1e-12` by default.
- `show_progress`: defines whether the progress informer must be displayed or not. `true` by default.
"""
macro evolution(args...)
    rules, loop = args[end-1:end]
    kwargs = Dict(a.args for a in args[begin:end-2])
    _evolution_block(rules, loop; kwargs...)
end
