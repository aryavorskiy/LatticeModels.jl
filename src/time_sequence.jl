import Base: ==

"""
    TimeSequence{ET}
Series of some data depending on time.
"""
struct TimeSequence{ET} <: AbstractDict{Float64, ET}
    times::Vector{Float64}
    snapshots::Vector{ET}
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
TimeSequence{ET}() where ET = TimeSequence{ET}([], [])
"""
    TimeSequence(value; timestamp=0)
Constructs a TimeSequence with one single snapshot. The timestamp is zero by default but can be overriden.
"""
TimeSequence(val; t::Real=0) = TimeSequence([t], [val])

const LatticeValueSequence = TimeSequence{LatticeValue}
const LatticeArraySequence = TimeSequence{LatticeArray}
const CurrentsSequence = TimeSequence{MaterializedCurrents}

"""
    time_domain(ts::TimeSequence)

Returns the timestamps of the snapshots.
"""
time_domain(tseq::TimeSequence) = tseq.times
snapshots(tseq::TimeSequence) = tseq.snapshots

Base.copy(tseq::TimeSequence) = TimeSequence(copy(tseq.times), [copy(s) for s in tseq.snapshots])
Base.empty(::TimeSequence{ET}) where ET = TimeSequence{ET}()
Base.length(tseq::TimeSequence) = length(time_domain(tseq))

function Base.show(io::IO, ::MIME"text/plain", tseq::TimeSequence{ET}) where ET
    print(io, "AbstractTimeSequence{$ET} with $(length(tseq)) records")
    td = time_domain(tseq)
    if length(td) ≥ 2
        print(io, "\nTimestamps in range $(td[1]) .. $(td[end])")
    elseif length(td) == 1
        print(io, "\nTimestamp: $(only(td))")
    end
end

==(ts1::TimeSequence, ts2::TimeSequence) =
    (ts1.times == ts2.times) && (ts1.snapshots == ts2.snapshots)

function Base.get(tseq::TimeSequence, t::Number, default)
    f = findfirst(==(t), tseq.times)
    f === nothing ? default : tseq.snapshots[f]
end

function Base.insert!(tseq::TimeSequence, t::Number, value)
    td = time_domain(tseq)
    imx = findlast(≤(t), td)
    imx === nothing && (imx = 0)
    if checkbounds(Bool, td, imx) && td[imx] == t
        tseq.snapshots[imx] = value
    else
        insert!(tseq.times, imx + 1, t)
        insert!(tseq.snapshots, imx + 1, value)
    end
    tseq
end

(tseq::TimeSequence)(t::Real) =
    tseq.snapshots[findmin(x -> abs(x - t), tseq.times)[2]]
function (tseq::TimeSequence)(tmin::Real, tmax::Real)
    inds = findall(time_domain(tseq)) do x
        tmin ≤ x ≤ tmax
    end
    TimeSequence(tseq.times[inds], tseq.snapshots[inds])
end

function Base.getindex(tseq::TimeSequence, args...)
    length(tseq) == 0 && error("Cannot index zero-length TimeSequence")
    TimeSequence(tseq.times, [tseq.snapshots[i][args...] for i in eachindex(tseq.snapshots)])
end

function Base.iterate(tseq::TimeSequence, state=(1, length(tseq)))
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    tseq.times[ind] => tseq.snapshots[ind], (ind + 1, len)
end

_internal(la::LatticeArray) = la.array
_internal(lv::LatticeValue) = lv.values
_internal(curr::MaterializedCurrents) = curr.currents
const LatticeType = Union{LatticeValue, LatticeArray, MaterializedCurrents}
function _axpby!(a, x::LatticeType, b, y::LatticeType)
    axpby!(a, _internal(x), b, _internal(y))
    y
end
_axpby!(a, x::Number, b, y::Number) = a * x + b * y
_axpby!(a, x, b, y) = axpby!(a, x, b, y)

"""
    differentiate!(ts::TimeSequence)

Differentiate the values stored in the `TimeSequence` object by time using the symmetric difference formula.
"""
function differentiate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot differentiate TimeSequence of length $(length(tseq))")
    td = time_domain(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i-1]
        tseq.snapshots[i-1] = _axpby!(1/dt, tseq.snapshots[i],
                -1/dt, tseq.snapshots[i - 1])
        td[i - 1] += dt / 2
    end
    pop!(td)
    pop!(tseq.snapshots)
    tseq
end
differentiate(tseq::TimeSequence) = differentiate!(copy(tseq))

"""
    integrate(ts::TimeSequence)

Integrates the values stored in the `TimeSequence` object over time using the trapezoidal rule.
"""
function integrate!(tseq::TimeSequence)
    length(tseq) < 2 && error("Cannot integrate TimeSequence of length $(length(tseq))")
    td = time_domain(tseq)
    for i in 2:length(tseq)
        dt = td[i] - td[i-1]
        tseq.snapshots[i-1] = _axpby!(1/2dt, tseq.snapshots[i], 1/2dt, tseq.snapshots[i-1])
    end
    last = pop!(tseq.snapshots)
    pushfirst!(tseq.snapshots, zero(last))
    for i in 2:length(tseq)
        tseq.snapshots[i] = _axpby!(1, tseq.snapshots[i-1], 1, tseq.snapshots[i])
    end
    tseq
end
integrate(tseq::TimeSequence) = integrate!(copy(tseq))
