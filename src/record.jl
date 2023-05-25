import Base: ==, diff

const STORABLE_TYPES = (LatticeArray, LatticeValue, MaterializedCurrents)
const StorableLatticeType = Union{STORABLE_TYPES...}
function _storable(ET::Type{<:StorableLatticeType})
    for NewET in STORABLE_TYPES
        ET <: NewET && return NewET
    end
    error("Non-storable type $ET cannot be a Record; expected one of $STORABLE_TYPES")
end

"""
    LatticeRecord{Eltype}

A struct providing the interface for storing information about
how some `LatticeValue`, `LatticeArray` or `MaterializedCurrents` depended on time.

Behaves like a vector of `(time, value)` tuples, supports time-based indexing (via call syntax, returns a stored record)
and site-based indexing (via bracket syntax, returns a vector or a new LatticeRecord, depending on the return type).
"""
struct LatticeRecord{ET<:StorableLatticeType} <: AbstractDict{Float64, ET}
    lattice::Lattice
    snapshots::Vector{Array}
    times::Vector{Float64}
    function LatticeRecord(vcs::AbstractVector, ts)
        !allequal(lattice.(vcs)) && throw(ArgumentError("all lattices must be equal"))
        new{_storable(eltype(vcs))}(lattice(first(vcs)), _internal.(vcs), ts)
    end
    LatticeRecord{ET}(l::Lattice, rcs, ts) where ET = new{_storable(ET)}(l, rcs, ts)
    LatticeRecord{ET}(l::Lattice) where ET = new{_storable(ET)}(l, [], [])
end

Base.first(lr::LatticeRecord{ET}) where ET = ET(lattice(lr), first(lr.snapshots))
Base.last(lr::LatticeRecord{ET}) where ET = ET(lattice(lr), last(lr.snapshots))
Base.length(lr::LatticeRecord) = length(lr.snapshots)

function Base.show(io::IO, ::MIME"text/plain", lr::LatticeRecord)
    print(io, "LatticeRecord with $(length(lr)) records")
    if length(lr.times) ≥ 2
        print(io, "\nTimestamps in range $(lr.times[1]) .. $(lr.times[end])")
    elseif length(lr.times) == 1
        println(io, "\nTimestamp: $(only(lr.times))")
    end
end

function ==(lr1::LatticeRecord{ET}, lr2::LatticeRecord{ET})  where {ET}
    (lr1.lattice == lr2.lattice) && (lr1.times == lr2.times) && (lr1.snapshots == lr2.snapshots)
end

const LatticeValueRecord = LatticeRecord{LatticeValue}
const LatticeArrayRecord = LatticeRecord{LatticeArray}
const CurrentsRecord = LatticeRecord{MaterializedCurrents}
_internal(la::LatticeArray) = la.array
_internal(lv::LatticeValue) = lv.values
_internal(curr::MaterializedCurrents) = curr.currents
init_record(slt::ET, t=0) where {ET<:StorableLatticeType} = LatticeRecord{ET}(lattice(slt), [_internal(slt)], Float64[t])

lattice(lr::LatticeRecord) = lr.lattice

"""
    time_domain(lr::LatticeRecord)

Returns the timestamps of the snapshots.
"""
time_domain(lr::LatticeRecord) = lr.times

function Base.insert!(lr::LatticeRecord{ET}, t::Number, value::ET) where ET
    check_lattice_match(lr, value)
    imx = findlast(≤(t), time_domain(lr))
    imx === nothing && (imx = 0)
    if checkbounds(Bool, lr.times, imx) && lr.times[imx] == t
        lr.snapshots[imx] = _internal(value)
    else
        insert!(lr.times, imx + 1, t)
        insert!(lr.snapshots, imx + 1, _internal(value))
    end
    lr
end

(lr::LatticeRecord{ET})(t::Real) where ET = ET(lr.lattice, lr.snapshots[findmin(x -> abs(x - t), lr.times)[2]])
function (lr::LatticeRecord{ET})(tmin::Real, tmax::Real) where ET
    inds = findall(time_domain(lr)) do x
        tmin ≤ x ≤ tmax
    end
    LatticeRecord{ET}(lattice(lr), lr.snapshots[inds], lr.times[inds])
end

function Base.getindex(lr::LatticeRecord{ET}, args...) where {ET}
    sample = first(lr)[args...]
    l = lattice(lr)
    if sample isa StorableLatticeType
        LatticeRecord([ET(l, rec)[args...] for rec in lr.snapshots], lr.times)
    else
        Dict(t => ET(lr.lattice, rec)[args...] for (t, rec) in zip(lr.times, lr.snapshots))
    end
end

function Base.iterate(lr::LatticeRecord{ET}, state=(1, length(lr))) where ET
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    lr.times[ind] => ET(lr.lattice, lr.snapshots[ind]), (ind + 1, len)
end

"""
    differentiate(lr::LatticeRecord)

Differentiate the values stored in the record by time using the symmetric difference formula.
"""
function differentiate(lr::LatticeRecord{ET}) where ET
    length(lr) < 2 && error("Cannot differentiate LatticeRecord of length $(length(lr))")
    td = time_domain(lr)
    LatticeRecord{ET}(lattice(lr),
        [@. (lr.snapshots[i + 1] - lr.snapshots[i])/(td[i + 1] - td[i])
        for i in 1:length(lr) - 1], (td[2:end] .+ td[1:end-1])./2)
end

"""
    integrate(lr::LatticeRecord)

Integrates the values stored in the record over time using the trapezoidal rule.
"""
function integrate(lr::LatticeRecord{ET}) where ET
    td = time_domain(lr)
    a = copy(first(lr.snapshots)) * ((td[2] - td[1]) / 2)
    out_snapshots = [zero(first(lr.snapshots))]
    for i in 2:length(lr)
        a .+= lr.snapshots[i] .* (td[i] - td[i-1])
        push!(out_snapshots, @. (a - lr.snapshots[i] / 2) * (td[i] - td[i-1]))
    end
    LatticeRecord{ET}(lattice(lr), out_snapshots, td)
end
