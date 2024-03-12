import Base: exp
using LinearAlgebra, ProgressMeter
import KrylovKit

"""
    SchroedingerSolver

Abstract type for solvers that can be used to evolve states in time according to the
Schroedinger equation.

See also concrete implementations: `CachedExp`, `KrylovKitExp`.

## Methods to implement
- `update_solver!(solver, hamiltonian, dt, force=false)`: Update the solver to evolve the
    states according to the given Hamiltonian and time step. If `force` is `true`, the
    solver should always update, even if the Hamiltonian and time step are seemingly the
    same as the previous ones.
- `step!(solver, state, cache)`: Evolve the given state in time using the solver. The `cache`
    argument is used to store intermediate results and can be `nothing` if the solver does
    not need it.
- `evolution_cache(solver, state)`: Return a cache object that can be used to store
    intermediate results for the given state. Returns `nothing` if the solver does not need
    a cache for the given state (this is the default implementation).
"""
abstract type SchroedingerSolver end
function update_solver! end
evolution_cache(solver::SchroedingerSolver, state::Ket) = evolution_cache(solver, state.data)
evolution_cache(solver::SchroedingerSolver, state::DataOperator) = evolution_cache(solver, state.data)
evolution_cache(::SchroedingerSolver, _) = nothing
step!(solver::SchroedingerSolver, state::Ket, cache) = step!(solver, state.data, cache)
step!(solver::SchroedingerSolver, state::DataOperator, cache) = step!(solver, state.data, cache)
step!(::SchroedingerSolver, ::Bra, ::Any) =
    throw(ArgumentError("Bra cannot be evolved in time; convert it to a Ket instead"))

const EvolutionHamType = Union{AbstractMatrix,DataOperator,AbstractTimeDependentOperator,Function}
const EvolutionStateType = Union{Ket,Bra,DataOperator,AbstractVecOrMat}
_data(arr::AbstractVecOrMat) = arr
_data(ket::Ket) = ket.data
_data(bra::Bra) = bra.data
_data(op::DataOperator) = op.data

"""
    CachedExp(ham[; threshold=1e-10, nonzero_tol=1e-14])

A `SchroedingerSolver` that finds the matrix exponential of the Hamiltonian and caches it.
The matrix exponential is computed using a scaling and squaring method, so this solver works
well with sparse or GPU arrays.

## Arguments
- `ham`: The Hamiltonian of the system. It can be an `Operator` or its matrix.
- `threshold`: The threshold for the error in the matrix exponential.
- `nonzero_tol`: The tolerance for dropping small elements in the matrix exponential if it is
    sparse.
"""
mutable struct CachedExp{MT,ET,KWT} <: SchroedingerSolver
    mat::MT
    matexp::ET
    factor::ComplexF64
    kw::KWT
end
CachedExp(ham; kw...) = CachedExp(zero(_data(ham)), one(complex(float(_data(ham)))), 0.0im, kw)
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
    solver.matexp = myexp!(solver.matexp, dmat, factor; solver.kw...)
    solver.factor = factor
end
function myexp!(P::AbstractMatrix, A::AbstractMatrix, factor; threshold=1e-10, nonzero_tol=1e-14)
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

"""
    KrylovKitExp(ham[; kw...])

A `SchroedingerSolver` that uses the `exponentiate` function from KrylovKit.jl to evolve the
wavefunction vectors. This solver is useful for large, sparse, time-dependent Hamiltonians.

## Arguments
- `ham`: The Hamiltonian of the system. It can be an `Operator` or its matrix.
- `kw...`: Keyword arguments to be passed to `exponentiate`.
"""
mutable struct KrylovKitExp{MT<:AbstractMatrix,KWT} <: SchroedingerSolver
    mat::MT
    factor::ComplexF64
    kw::KWT
end
KrylovKitExp(ham; kw...) = KrylovKitExp(zero(_data(ham)), 0.0im, kw)
function update_solver!(solver::KrylovKitExp, mat, dt, _force)
    solver.factor = -im * dt
    solver.mat = _data(mat)
end
function step!(solver::KrylovKitExp, state::AbstractVector, _cache)
    newstate, info = KrylovKit.exponentiate(solver.mat, solver.factor, state; solver.kw...)
    info.converged == 0 && throw(ArgumentError("`exponentiate` did not converge"))
    copyto!(state, newstate)
end

"""
    Evolution([solver, ]hamiltonian, states...; timedomain, namedstates...)

Create an `Evolution` object that can be used to evolve states in time according to the
Schroedinger equation.

# Arguments
- `solver`: A `SchroedingerSolver` object that will be used to evolve the states. If omitted,
    a `CachedExp` solver will be created.
- `hamiltonian`: The Hamiltonian of the system. It can be a matrix, a time-dependent operator
    or a function that returns the Hamiltonian at a given time.
- `states` and `namedstates`: The states to be evolved. They can be `Ket` wavefunctions or
    `DataOperator` density matrices.
- `timedomain`: The time domain to be used for the evolution. If omitted, the non-iterable `Evolution`
    object will be returned, and you will be able to call it with
    the time domain later.

See `SchroedingerSolver` for more information about solvers.

!!! warning
    Please note that the `Evolution` object is a **stateful** iterator. This means that it
    keeps track of the current time and the states as they evolve. You can perform evolution
    multiple times, but the timeline will be kept and the states will be updated in place.

    Also do not edit the states in place, as this will affect the evolution. If you need to
    modify the states or save them, make a copy of them first.
"""
struct Evolution{SolverT,HamT,NamedTupleT}
    solver::SolverT
    hamiltonian::HamT
    states::NamedTupleT
    time::Base.RefValue{Float64}
    function Evolution(solver::SolverT, hamiltonian::HamT, states::Union{Tuple,NamedTuple}) where
            {SolverT<:SchroedingerSolver,HamT}
        newstates = map(states) do state
            state isa EvolutionStateType || throw(ArgumentError("invalid state type: $(typeof(state))"))
            cache = evolution_cache(solver, state)
            return copy(state), cache
        end
        return new{SolverT,HamT,typeof(newstates)}(solver, hamiltonian, newstates, Ref(0.0))
    end
end
function Evolution(hamiltonian::EvolutionHamType, states::EvolutionStateType...; namedstates...)
    solver = CachedExp(eval_hamiltonian(hamiltonian, 0))
    return Evolution(solver, hamiltonian, states...; namedstates...)
end
function Evolution(solver::SchroedingerSolver, hamiltonian::EvolutionHamType, states::EvolutionStateType...;
        timedomain=nothing, namedstates...)
    final_states = if isempty(states)
        isempty(namedstates) && throw(ArgumentError("No states provided"))
        NamedTuple(namedstates)
    elseif isempty(namedstates)
        states
    else
        throw(ArgumentError("Do not use named and unnamed states together"))
    end
    evol = Evolution(solver, hamiltonian, final_states)
    if timedomain === nothing
        return evol
    else
        return evol(timedomain)
    end
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
Base.length(iter::EvolutionIterator) = length(iter.times)
function Base.iterate(iter::EvolutionIterator)
    p = Progress(length(iter), dt=0.3, desc="Unitary evolution... ", showspeed=true,
        barglyphs=BarGlyphs("[=> ]"))
    return iterate(iter, (1, p))
end
function Base.iterate(iter::EvolutionIterator, st)
    i, p = st
    i > length(iter.times) && return nothing
    dt = i > 1 ? iter.times[i] - iter.times[i-1] : iter.times[i]
    t = iter.evol.time[]
    ham = step!(iter.evol, dt)
    states = map(first, iter.evol.states)
    next!(p)
    return EvolutionTimestamp(ham, states, Float64(t)), (i+1, p)
end

struct EvolutionTimestamp{HamT, NamedTupleT}
    H::HamT
    states::NamedTupleT
    t::Float64
end
Base.iterate(ts::EvolutionTimestamp, i::Int=1) =
    i > length(ts.states) ? (ts.H, Val(:time)) : (ts.states[i], i+1)
Base.iterate(ts::EvolutionTimestamp, ::Val{:time}) = (ts.t, Val(:end))
Base.iterate(::EvolutionTimestamp, ::Val{:end}) = nothing

Base.getindex(ts::EvolutionTimestamp, sym::Symbol) = ts.states[sym]
Base.getindex(ts::EvolutionTimestamp, i::Int) = ts.states[i]

@inline function Base.getproperty(ts::EvolutionTimestamp, sym::Symbol)
    if sym in fieldnames(EvolutionTimestamp)
        return getfield(ts, sym)
    elseif sym == :state
        return only(getfield(ts, :states))
    else
        return getfield(ts, :states)[sym]
    end
end
