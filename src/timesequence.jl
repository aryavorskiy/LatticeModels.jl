"""
    TimeSequence{ET}

A time-ordered sequence of values.
"""
struct TimeSequence{ET} <: AbstractDict{Float64, ET}
    times::Vector{Float64}
    values::Vector{ET}
    function TimeSequence{ET}(ts, vs) where ET
        length(ts) != length(vs) &&
            throw(ArgumentError("Keys/values length mismatch:\n$(length(ts)) timestamps, $(length(vs)) snapshots"))
        new{ET}(ts, vs)
    end
end
"""
    TimeSequence(times, values)

Constructs a `TimeSequence` with the given `times` and `values`.
"""
TimeSequence(ts, vs::AbstractVector) = TimeSequence{eltype(vs)}(ts, vs)

"""
    TimeSequence{ET}()

Constructs an empty `TimeSequence` with eltype `ET`.
"""
TimeSequence{ET}() where ET = TimeSequence{ET}(Float64[], ET[])

"""
    TimeSequence(f, ev, times)
    TimeSequence(f, ev_iter)

Constructs a TimeSequence by iterating the evolution iterator `ev` and applying the function `f` to each moment.

## Arguments
- `f`: A function that takes a moment and returns a value. The function is applied to each moment in the evolution.
- `ev`: An `Evolution` object.
- `times`: A range of times to evaluate the function at.
- `ev_iter`: An evoltuion iterator that yields moments. Think of it as `ev(times)`.
"""
TimeSequence(f::Function, evol::EvolutionIterator) =
    TimeSequence(evol.times, [f(moment) for moment in evol])
TimeSequence(f::Function, evol::Evolution, times) = TimeSequence(f, evol(times))

"""
    timestamps(tseq::TimeSequence)

Returns the timestamps of the `TimeSequence`.
"""
timestamps(tseq::TimeSequence) = tseq.times
"""
    timerange(tseq::TimeSequence)

Returns the range of the timestamps of the `TimeSequence`.
"""
timerange(tseq::TimeSequence) = first(tseq.times)..last(tseq.times)
Base.values(tseq::TimeSequence) = tseq.values

Base.copy(tseq::TimeSequence) = TimeSequence(copy(tseq.times), [copy(s) for s in tseq.values])
Base.length(tseq::TimeSequence) = length(timestamps(tseq))
Base.map(f, tseq::TimeSequence) = TimeSequence(timestamps(tseq), map(f, tseq.values))

function Base.show(io::IO, mime::MIME"text/plain", tseq::TimeSequence{ET}) where ET
    print(io, "TimeSequence{$ET} with ", fmtnum(tseq, "entr", "y", "ies"))
    requires_compact(io) && return
    length(tseq) ≥ 2 && print(io, "\nTimestamps in range $(timerange(tseq)):")
    maxlen = get(io, :maxlines, 10)
    io = IOContext(io, :compact => true)
    for i in 1:min(length(tseq), maxlen)
        print(io, "\n")
        if i == maxlen < length(tseq)
            print(io, "  ⋮")
        else
            print(io, "  $(round(tseq.times[i], digits=5)) => ")
            show(io, mime, tseq.values[i])
        end
    end
end

Base.:(==)(ts1::TimeSequence, ts2::TimeSequence) =
    (ts1.times == ts2.times) && (ts1.values == ts2.values)

Base.empty(::TimeSequence{ET}) where ET = TimeSequence{ET}()
Base.empty(::TimeSequence, ::Type, ::VT) where VT =
    TimeSequence(VT)()
function Base.get(tseq::TimeSequence, t::Real, default)
    f = findfirst(≈(t, atol=√eps()), tseq.times)
    f === nothing ? default : tseq.values[f]
end
function Base.setindex!(tseq::TimeSequence, val, t::Real)
    i = findfirst(≈(t, atol=√eps()), tseq.times)
    if i === nothing
        i = searchsortedfirst(tseq.times, t)
        insert!(tseq.times, i, t)
        insert!(tseq.values, i, val)
    else
        tseq.values[i] = val
    end
    return val
end
function Base.delete!(tseq::TimeSequence, t::Number)
    f = findfirst(≈(t, atol=√eps()), tseq.times)
    f === nothing && return
    deleteat!(tseq.times, f)
    deleteat!(tseq.values, f)
    return tseq
end

Base.getindex(tseq::TimeSequence, t::Real) = tseq[t = t]
@inline _index_inner(val, ::Tuple{}, ::Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, @NamedTuple{}}) = val
@inline _index_inner(val, args::Tuple, kw::Base.Iterators.Pairs) = getindex(val, args...; kw...)
function Base.getindex(tseq::TimeSequence, args...; t=-Inf..Inf, kw...)
    if t isa Number
        val = get(tseq, t, Base.secret_table_token)
        if val == Base.secret_table_token
            throw(KeyError(t))
        else
            return _index_inner(val, args, kw)
        end
    else
        mask = map(in(t), tseq.times)
        return TimeSequence(tseq.times[mask],
        [_index_inner(tseq.values[i], args, kw) for i in eachindex(tseq.values) if mask[i]])
    end
end
function Base.delete!(tseq::TimeSequence; t)
    f = findall(in(t), tseq.times)
    deleteat!(tseq.times, f)
    deleteat!(tseq.values, f)
    return tseq
end

function Base.iterate(tseq::TimeSequence, state=(1, length(tseq)))
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    return tseq.times[ind] => tseq.values[ind], (ind + 1, len)
end

_underlying_array(op::DataOperator) = op.data
_underlying_array(lv::LatticeValue) = lv.values
_underlying_array(curr::Currents) = curr.currents
_underlying_array(arr::AbstractArray) = arr
_underlying_array(::Any) = nothing
_underlying_array(::StaticArray) = nothing
function _axpby!(a, x::T, b, y::T) where T
    if _underlying_array(x) === nothing
        return a * x + b * y
    else
        axpby!(a, _underlying_array(x), b, _underlying_array(y))
        return y
    end
end

"""
    differentiate!(tseq::TimeSequence)

Differentiate the values stored in `tseq` by time using the symmetric difference formula.
The new values are written into `tseq`.

## Example
```jldoctest
julia> using LatticeModels

julia> tseq = TimeSequence(0:0.1:10, 0:0.1:10)  # f(t) = t
TimeSequence{Float64} with 101 entry
Timestamps in range 0.0 .. 10.0:
  0.0 => 0.0
  0.1 => 0.1
  0.2 => 0.2
  0.3 => 0.3
  0.4 => 0.4
  0.5 => 0.5
  0.6 => 0.6
  0.7 => 0.7
  0.8 => 0.8
  ⋮

julia> differentiate!(tseq)                     # f'(t) = 1
TimeSequence{Float64} with 100 entries
Timestamps in range 0.05 .. 9.95:
  0.05 => 1.0
  0.15 => 1.0
  0.25 => 1.0
  0.35 => 1.0
  0.45 => 1.0
  0.55 => 1.0
  0.65 => 1.0
  0.75 => 1.0
  0.85 => 1.0
  ⋮
```
"""
function differentiate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot differentiate TimeSequence of length $(length(tseq))")
    td = timestamps(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i - 1]
        tseq.values[i - 1] = _axpby!(1/dt, tseq.values[i],
                -1/dt, tseq.values[i - 1])
        td[i - 1] += dt / 2
    end
    pop!(td)
    pop!(tseq.values)
    tseq
end
"""
    differentiate(tseq::TimeSequence)

Differentiate the values stored in `tseq` and create a copy; see [`differentiate!`](@ref).
"""
differentiate(tseq::TimeSequence) = differentiate!(copy(tseq))

"""
    integrate!(tseq::TimeSequence)

Integrate the values stored in `tseq` over time using the trapezoidal rule.
The new values are written into `tseq`. The first value is set to zero.

## Example
```jldoctest
julia> using LatticeModels

julia> tseq = TimeSequence(0:0.1:10, 0:0.1:10)  # f(t) = t
TimeSequence{Float64} with 101 entry
Timestamps in range 0.0 .. 10.0:
  0.0 => 0.0
  0.1 => 0.1
  0.2 => 0.2
  0.3 => 0.3
  0.4 => 0.4
  0.5 => 0.5
  0.6 => 0.6
  0.7 => 0.7
  0.8 => 0.8
  ⋮

julia> integrate!(tseq)                         # F(t) = t^2 / 2
TimeSequence{Float64} with 101 entry
Timestamps in range 0.0 .. 10.0:
  0.0 => 0.0
  0.1 => 0.005
  0.2 => 0.02
  0.3 => 0.045
  0.4 => 0.08
  0.5 => 0.125
  0.6 => 0.18
  0.7 => 0.245
  0.8 => 0.32
  ⋮
```
"""
function integrate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot integrate TimeSequence of length $(length(tseq))")
    td = timestamps(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i-1]
        tseq.values[i-1] = _axpby!(dt/2, tseq.values[i], dt/2, tseq.values[i-1])
    end
    last = pop!(tseq.values)
    pushfirst!(tseq.values, zero(last))
    for i in 2:length(tseq)
        tseq.values[i] = _axpby!(1, tseq.values[i-1], 1, tseq.values[i])
    end
    tseq
end
"""
    integrate(tseq::TimeSequence)

Integrate the values stored in `tseq` and create a copy; see [`integrate!`](@ref).
"""
integrate(tseq::TimeSequence) = integrate!(copy(tseq))
