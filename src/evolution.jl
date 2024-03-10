import Base: exp
using LinearAlgebra, ProgressMeter
using KrylovKit

abstract type SchrodingerSolver end
function update_solver! end
evolution_cache(solver::SchrodingerSolver, state::Ket) = evolution_cache(solver, state.data)
evolution_cache(solver::SchrodingerSolver, state::DataOperator) = evolution_cache(solver, state.data)
evolution_cache(::SchrodingerSolver, _) = nothing
step!(solver::SchrodingerSolver, state::Ket, cache) = step!(solver, state.data, cache)
step!(solver::SchrodingerSolver, state::DataOperator, cache) = step!(solver, state.data, cache)
step!(::SchrodingerSolver, ::Bra, ::Any) =
    throw(ArgumentError("Bra cannot be evolved in time; convert it to a Ket instead"))

const EvolutionStateType = Union{Ket,Bra,DataOperator,AbstractVecOrMat}
_data(arr::AbstractVecOrMat) = arr
_data(ket::Ket) = ket.data
_data(bra::Bra) = bra.data
_data(op::DataOperator) = op.data
mutable struct CachedExp{MT,ET} <: SchrodingerSolver
    mat::MT
    matexp::ET
    factor::ComplexF64
end
CachedExp(mat) = CachedExp(zero(_data(mat)), one(complex(float(_data(mat)))), 0.0im)
function step!(solver::CachedExp, state::AbstractVector, cache)
    copyto!(cache, state)
    mul!(state, solver.matexp, cache)
end
function step!(solver::CachedExp, state::AbstractMatrix, cache)
    c1, c2 = cache
    copyto!(c1, state)
    mul!(c2, solver.matexp, c1)
    mul!(state, c2, solver.matexp')
end
evolution_cache(::CachedExp, state::AbstractVector) = similar(state)
evolution_cache(::CachedExp, state::AbstractMatrix) = (similar(state), similar(state))

function update_solver!(solver::CachedExp, mat, dt, force=false)
    factor = -im * dt
    dmat = _data(mat)
    if !force
        factor ≈ solver.factor && (solver.mat === dmat || solver.mat == dmat) && return
    end
    solver.mat = dmat
    solver.matexp = myexp!(solver.matexp, dmat, factor)
    solver.factor = factor
end
function myexp!(P::AbstractMatrix, A::AbstractMatrix, factor; threshold=1e-8, nonzero_tol=1e-14)
    mat_norm = norm(A, Inf)
    (iszero(mat_norm) || iszero(factor)) && return one(A)
    scaling_factor = nextpow(2, mat_norm * abs(factor))
    A /= scaling_factor
    delta = one(mat_norm)
    @check_size A :square

    P .= one(A)
    next_term = copy(P)
    nt_buffer = similar(next_term)
    n = 1
    while delta > threshold / scaling_factor
        if issparse(A)
            next_term = A * next_term
            next_term .*= factor / n
            droptol!(next_term, nonzero_tol)
        else
            mul!(nt_buffer, A, next_term, factor / n, 0)
            next_term, nt_buffer = nt_buffer, next_term
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
            P, nt_buffer = nt_buffer, P
        end
    end
    return P
end

struct Evolution{SolverT,HamT,NamedTupleT}
    solver::SolverT
    hamiltonian::HamT
    states::NamedTupleT
    time::Base.RefValue{Float64}
    function Evolution(solver::SolverT, hamiltonian::HamT, states::NamedTuple) where
            {SolverT<:SchrodingerSolver,HamT}
        newstates = map(states) do state
            state isa EvolutionStateType || throw(ArgumentError("invalid state type: $(typeof(state))"))
            cache = evolution_cache(solver, state)
            return copy(state), cache
        end
        return new{SolverT,HamT,typeof(newstates)}(solver, hamiltonian, newstates, Ref(0.0))
    end
end
function Evolution(solver::SchrodingerSolver, hamiltonian; kw...)
    return Evolution(solver, hamiltonian, NamedTuple(kw))
end
function Evolution(hamiltonian; kw...)
    solver = CachedExp(eval_hamiltonian(hamiltonian, 0))
    return Evolution(solver, hamiltonian, NamedTuple(kw))
end

eval_hamiltonian(hamiltonian, _) = hamiltonian
eval_hamiltonian(hamiltonian::Function, t) = hamiltonian(t)
function eval_hamiltonian(ham::QuantumOpticsBase.AbstractTimeDependentOperator, t)
    QuantumOpticsBase.set_time!(ham, t)
    return ham
end
function step!(evol::Evolution, dt)
    t = evol.time[]
    H = eval_hamiltonian(evol.hamiltonian, t)
    abs(dt) < 1e-15 && return H
    update_solver!(evol.solver, H, dt, H isa QuantumOpticsBase.AbstractTimeDependentOperator)
    for (state, cache) in values(evol.states)
        step!(evol.solver, state, cache)
    end
    evol.time[] += dt
    return H
end

(evol::Evolution)(ts::AbstractVector{<:Real}) = EvolutionIterator(evol, ts)

struct EvolutionIterator{EvolutionT,TimesT}
    evol::EvolutionT
    times::TimesT
end
@inline function Base.iterate(iter::EvolutionIterator, i=1)
    i > length(iter.times) && return nothing
    dt = i > 1 ? iter.times[i] - iter.times[i-1] : iter.times[i]
    t = iter.evol.time[]
    ham = step!(iter.evol, dt)
    states = map(first, iter.evol.states)
    return states, i+1
    return EvolutionTimestamp(ham, states, Float64(t)), i+1
end

struct EvolutionTimestamp{HamT, NamedTupleT}
    H::HamT
    states::NamedTupleT
    t::Float64
end
Base.iterate(ts::EvolutionTimestamp) = ts.H, Val(:time)
Base.iterate(ts::EvolutionTimestamp, ::Val{:time}) = ts.t, 1
Base.iterate(ts::EvolutionTimestamp, i::Int) =
    i > length(ts.states) ? nothing : ts.states[i], i+1

Base.getindex(ts::EvolutionTimestamp, sym::Symbol) = ts.states[sym]
Base.getindex(ts::EvolutionTimestamp, i::Int) = ts.states[i]
