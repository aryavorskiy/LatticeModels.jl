import Base: ==

"""
    time_domain(ts::TimeSequence)

Returns the timestamps of the snapshots.
"""
function time_domain end
function _internals end

struct TimeSequence{ET} <: AbstractDict{Float64, ET}
    times::Vector{Float64}
    snapshots::Vector{ET}
    function TimeSequence{ET}(times, snapshots) where ET
        length(times) != length(snapshots) &&
            error("Keys/values length mismatch:\n$(length(times)) timestamps, $(length(snapshots)) snapshots")
        new{ET}(times, snapshots)
    end
end
TimeSequence{ET}() where ET = TimeSequence{ET}([], [])
TimeSequence(ts, vcs::AbstractVector) = TimeSequence{eltype(vcs)}(ts, vcs)
TimeSequence(val; t::Real=0) = TimeSequence([t], [val])

const LatticeValueSequence = TimeSequence{LatticeValue}
const LatticeArraySequence = TimeSequence{LatticeArray}
const CurrentsSequence = TimeSequence{MaterializedCurrents}

time_domain(lr::TimeSequence) = lr.times
snapshots(lr::TimeSequence) = lr.snapshots

Base.copy(ts::TimeSequence) = TimeSequence(copy(ts.times), [copy(s) for s in snapshots(ts)])
Base.empty(::TimeSequence{ET}) where ET = TimeSequence{ET}()
Base.length(lr::TimeSequence) = length(time_domain(lr))

function Base.show(io::IO, ::MIME"text/plain", lr::TimeSequence{ET}) where ET
    print(io, "AbstractTimeSequence{$ET} with $(length(lr)) records")
    times = time_domain(lr)
    if length(times) ≥ 2
        print(io, "\nTimestamps in range $(times[1]) .. $(times[end])")
    elseif length(times) == 1
        print(io, "\nTimestamp: $(only(times))")
    end
end

==(lr1::TimeSequence{ET}, lr2::TimeSequence{ET})  where {ET} =
    (lr1.times == lr2.times) && (lr1.snapshots == lr2.snapshots)

function Base.insert!(lr::TimeSequence{ET}, t::Number, value::ET) where ET
    times = time_domain(lr)
    imx = findlast(≤(t), times)
    imx === nothing && (imx = 0)
    if checkbounds(Bool, times, imx) && times[imx] == t
        lr.snapshots[imx] = value
    else
        insert!(lr.times, imx + 1, t)
        insert!(lr.snapshots, imx + 1, value)
    end
    lr
end

(lr::TimeSequence{ET})(t::Real) where ET =
    lr.snapshots[findmin(x -> abs(x - t), lr.times)[2]]
function (lr::TimeSequence{ET})(tmin::Real, tmax::Real) where ET
    inds = findall(time_domain(lr)) do x
        tmin ≤ x ≤ tmax
    end
    TimeSequence(lr.times[inds], lr.snapshots[inds])
end

function Base.getindex(lr::TimeSequence, args...)
    length(lr) == 0 && error("Cannot index zero-length TimeSequence")
    TimeSequence(lr.times, [lr.snapshots[i][args...] for i in eachindex(lr.snapshots)])
end

function Base.iterate(lr::TimeSequence{ET}, state=(1, length(lr))) where ET
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    lr.times[ind] => lr.snapshots[ind], (ind + 1, len)
end

_internal(la::LatticeArray) = la.array
_internal(lv::LatticeValue) = lv.values
_internal(curr::MaterializedCurrents) = curr.currents
const LatticeType = Union{LatticeValue, LatticeArray, MaterializedCurrents}
function _axpby!(a, x::LatticeType, b, y::LatticeType)
    axpby!(a, _internal(x), b, _internal(y))
    y
end
_axpby!(a, x::AbstractArray, b, y::AbstractArray) = axpby!(a, x, b, y)
_axpby!(a, x::Number, b, y::Number) = a * x + b * y

"""
    differentiate!(ts::TimeSequence)

Differentiate the values stored in the record by time using the symmetric difference formula.
"""
function differentiate!(ts::TimeSequence)
    length(ts) < 2 && error("Cannot differentiate TimeSequence of length $(length(ts))")
    td = time_domain(ts)
    for i in 2:length(ts)
        dt = td[i] - td[i-1]
        ts.snapshots[i-1] = _axpby!(1/dt, ts.snapshots[i],
                -1/dt, ts.snapshots[i - 1])
        td[i - 1] += dt / 2
    end
    pop!(td)
    pop!(ts.snapshots)
    ts
end
differentiate(ts::TimeSequence) = differentiate!(copy(ts))

"""
    integrate(ts::TimeSequence)

Integrates the values stored in the record over time using the trapezoidal rule.
"""
function integrate!(ts::TimeSequence)
    length(ts) < 2 && error("Cannot integrate TimeSequence of length $(length(ts))")
    td = time_domain(ts)
    for i in 2:length(ts)
        dt = td[i] - td[i-1]
        ts.snapshots[i-1] = _axpby!(1/2dt, ts.snapshots[i], 1/2dt, ts.snapshots[i-1])
    end
    last = pop!(ts.snapshots)
    pushfirst!(ts.snapshots, zero(last))
    for i in 2:length(ts)
        ts.snapshots[i] = _axpby!(1, ts.snapshots[i-1], 1, ts.snapshots[i])
    end
    ts
end
integrate(ts::TimeSequence) = integrate!(copy(ts))
