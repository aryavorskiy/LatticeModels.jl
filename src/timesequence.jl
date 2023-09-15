"""
    TimeSequence{ET}
Series of some data depending on time.
"""
struct TimeSequence{ET} <: AbstractDict{Float64, ET}
    times::Vector{Float64}
    values::Vector{ET}
    function TimeSequence{ET}(ts, vs) where ET
        length(ts) != length(vs) &&
            error("Keys/values length mismatch:\n$(length(ts)) timestamps, $(length(vs)) snapshots")
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
    TimeSequence(value; t=0)
Constructs a TimeSequence with one single snapshot. The timestamp is zero by default but can be over
riden.
"""
TimeSequence(t::Real, val) = TimeSequence(Float64[t], [val])

"""
    timestamps(ts::TimeSequence)

Returns the timestamps of the snapshots.
"""
timestamps(tseq::TimeSequence) = tseq.times
timerange(tseq::TimeSequence) = first(tseq.times)..last(tseq.times)
snapshots(tseq::TimeSequence) = tseq.values

Base.copy(tseq::TimeSequence) = TimeSequence(copy(tseq.times), [copy(s) for s in tseq.values])
Base.empty(::TimeSequence{ET}) where ET = TimeSequence{ET}()
Base.length(tseq::TimeSequence) = length(timestamps(tseq))

function Base.show(io::IO, ::MIME"text/plain", tseq::TimeSequence{ET}) where ET
    print(io, "AbstractTimeSequence{$ET} with $(length(tseq)) records")
    td = timestamps(tseq)
    if length(td) ≥ 2
        print(io, "\nTimestamps in range $(timerange(tseq))")
    elseif length(td) == 1
        print(io, "\nTimestamp: $(only(td))")
    end
end

Base.:(==)(ts1::TimeSequence, ts2::TimeSequence) =
    (ts1.times == ts2.times) && (ts1.values == ts2.values)

Base.empty(::TimeSequence, ::Type, ::VT) where VT =
    TimeSequence(VT)()
function Base.get(tseq::TimeSequence, t::Real, default)
    f = findfirst(≈(t, atol=√eps()), tseq.times)
    f === nothing ? default : tseq.values[f]
end
get_inner(val, ::Nothing) = val
get_inner(val, tup::Tuple) = val[tup...]
get_inner(val, ind) = val[ind]
function Base.getindex(ts::TimeSequence, domain; inner=nothing)
    mask = map(in(domain), ts.times)
    return TimeSequence(ts.times[mask],
        [get_inner(ts.values[i], inner) for i in eachindex(ts.values) if mask[i]])
end
function Base.getindex(ts::TimeSequence, t::Number; inner=nothing)
    val = get(ts, t, Base.secret_table_token)
    if val == Base.secret_table_token
        throw(KeyError(t))
    else
        return get_inner(val, inner)
    end
end
Base.getindex(ts::TimeSequence; inner) = ts[timerange(ts), inner=inner]
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
function Base.delete!(tseq::TimeSequence, t)
    f = findfirst(≈(t, atol=√eps()), tseq.times)
    f === nothing && return
    deleteat!(tseq.times, f)
    deleteat!(tseq.values, f)
end

function Base.iterate(tseq::TimeSequence, state=(1, length(tseq)))
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    tseq.times[ind] => tseq.values[ind], (ind + 1, len)
end

_internal(op::DataOperator) = op.data
_internal(lv::LatticeValue) = lv.values
_internal(curr::MaterializedCurrents) = curr.currents
# const LatticeType = Union{LatticeValue, LatticeArray, MaterializedCurrents}
function _axpby!(a, x, b, y)
    axpby!(a, _internal(x), b, _internal(y))
    y
end
_axpby!(a, x::Number, b, y::Number) = a * x + b * y
# _axpby!(a, x, b, y) = axpby!(a, x, b, y)

"""
    differentiate!(ts::TimeSequence)

Differentiate the values stored in the `TimeSequence` object by time using the symmetric difference formula.
"""
function differentiate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot differentiate TimeSequence of length $(length(tseq))")
    td = timestamps(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i-1]
        tseq.values[i-1] = _axpby!(1/dt, tseq.values[i],
                -1/dt, tseq.values[i - 1])
        td[i - 1] += dt / 2
    end
    pop!(td)
    pop!(tseq.values)
    tseq
end
differentiate(tseq::TimeSequence) = differentiate!(copy(tseq))

"""
    integrate(ts::TimeSequence)

Integrates the values stored in the `TimeSequence` object over time using the trapezoidal rule.
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
