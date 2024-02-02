"""
    TimeSequence{ET}
Series of some data depending on time.
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
Constructs a TimeSequence with given timestamps and values.
"""
TimeSequence(ts, vs::AbstractVector) = TimeSequence{eltype(vs)}(ts, vs)
TimeSequence{ET}() where ET = TimeSequence{ET}(Float64[], ET[])
"""
    TimeSequence(value[; t=0])
Constructs a TimeSequence with one single snapshot.
The timestamp is zero by default but can be over
riden.
"""
TimeSequence(t::Real, val) = TimeSequence(Float64[t], [val])

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

function Base.show(io::IO, ::MIME"text/plain", tseq::TimeSequence{ET}) where ET
    print(io, "TimeSequence{$ET} with $(length(tseq)) records")
    td = timestamps(tseq)
    if length(td) ≥ 2
        print(io, "\nTimestamps in range $(timerange(tseq))")
    elseif length(td) == 1
        print(io, "\nTimestamp: $(only(td))")
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
@inline _index_inner(val, ::Tuple{}) = val
@inline _index_inner(val, args::Tuple) = val[args...]
function Base.getindex(tseq::TimeSequence, args...; t)
    if t isa Number
        val = get(tseq, t, Base.secret_table_token)
        if val == Base.secret_table_token
            throw(KeyError(t))
        else
            return _index_inner(val, args)
        end
    else
        mask = map(in(t), tseq.times)
        return TimeSequence(tseq.times[mask],
        [_index_inner(tseq.values[i], args) for i in eachindex(tseq.values) if mask[i]])
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
The new values are written into `tseq`.
"""
function integrate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot integrate TimeSequence of length $(length(tseq))")
    td = timestamps(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i-1]
        tseq.values[i-1] = _axpby!(1/2dt, tseq.values[i], 1/2dt, tseq.values[i-1])
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

@recipe function f(tseq::TimeSequence)
    tseq.times, tseq.values
end

mutable struct TimeSequenceContainer
    seq::Nullable{TimeSequence}
    TimeSequenceContainer() = new(nothing)
end

function Base.setindex!(tsc::TimeSequenceContainer, val::T, t::Real) where T
    if tsc.seq === nothing
        tsc.seq = TimeSequence{T}([t], [val])
    else
        tsc.seq[t] = val
    end
    return val
end
