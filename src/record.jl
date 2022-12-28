import Base: ==, length, first, last, insert!, diff

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
struct LatticeRecord{ET<:StorableLatticeType} <: Function
    lattice::Lattice
    records::Vector{Array}
    times::Vector{Float64}
    function LatticeRecord(vcs::Vector{ET}, ts) where {ET}
        !allequal(lattice.(vcs)) && throw(ArgumentError("all lattices must be equal"))
        new{_storable(ET)}(lattice(first(vcs)), _internal.(vcs), ts)
    end
    LatticeRecord{ET}(l::Lattice, rcs, ts) where ET = new{_storable(ET)}(l, rcs, ts)
    LatticeRecord{ET}(l::Lattice) where ET = new{_storable(ET)}(l, [], [])
end

first(lvr::LatticeRecord{ET}) where ET = ET(lattice(lvr), first(lvr.records))
last(lvr::LatticeRecord{ET}) where ET = ET(lattice(lvr), last(lvr.records))
length(lvr::LatticeRecord) = length(lvr.records)

function show(io::IO, ::MIME"text/plain", lr::LatticeRecord)
    print(io, "LatticeRecord with $(length(lr)) records")
    if length(lr.times) ≥ 2
        print(io, "\nTimestamps in range $(lr.times[1]) .. $(lr.times[end])")
    elseif length(lr.times) == 1
        println(io, "\nTimestamp: $(only(lr.times))")
    end
end

function ==(lr1::LatticeRecord{ET}, lr2::LatticeRecord{ET})  where {ET}
    (lr1.lattice == lr2.lattice) && (lr1.times == lr2.times) && (lr1.records == lr2.records)
end

const LatticeValueRecord = LatticeRecord{LatticeValue}
const LatticeArrayRecord = LatticeRecord{LatticeArray}
const CurrentsRecord = LatticeRecord{MaterializedCurrents}
_internal(la::LatticeArray) = la.array
_internal(lv::LatticeValue) = lv.values
_internal(curr::MaterializedCurrents) = curr.currents
init_record(slt::ET) where {ET<:StorableLatticeType} = LatticeRecord{ET}(lattice(slt), [_internal(slt)], [0.])

lattice(lvr::LatticeRecord) = lvr.lattice

"""
    time_domain(lr::LatticeRecord)

Returns the timestamps of the records.
"""
time_domain(lr::LatticeRecord) = lr.times

function Base.insert!(lvr::LatticeRecord{ET}, t::Number, value::ET) where ET
    check_lattice_match(lvr, value)
    imx = findlast(≤(t), time_domain(lvr))
    imx === nothing && (imx = 0)
    if checkbounds(Bool, lvr.times, imx) && lvr.times[imx] == t
        lvr.records[imx] = _internal(value)
    else
        insert!(lvr.times, imx + 1, t)
        insert!(lvr.records, imx + 1, _internal(value))
    end
    lvr
end

(lr::LatticeRecord{ET})(t::Real) where ET = ET(lr.lattice, lr.records[findmin(x -> abs(x - t), lr.times)[2]])
function (lr::LatticeRecord{ET})(tmin::Real, tmax::Real) where ET
    inds = findall(time_domain(lr)) do x
        tmin ≤ x ≤ tmax
    end
    LatticeRecord{ET}(lattice(lr), lr.records[inds], lr.times[inds])
end
function Base.getindex(lr::LatticeRecord{ET}, args...) where {ET}
    sample = first(lr)[args...]
    l = lattice(lr)
    if sample isa StorableLatticeType
        LatticeRecord([ET(l, rec)[args...] for rec in lr.records], lr.times)
    else
        [ET(l, rec)[args...] for rec in lr.records]
    end
end

function iterate(lvr::LatticeRecord{ET}, state=(1, length(lattice(lvr)))) where ET
    ind, len = state
    1 ≤ ind ≤ len || return nothing
    (lvr.times[ind], ET(lvr.lattice, lvr.records[ind])), (ind + 1, len)
end

"""
    diff(lattice_record)

Differentiate the values stored in the record by time using the symmetric difference formula.
"""
function Base.diff(lvr::LatticeRecord{ET}) where ET
    length(lvr) < 2 && error("Cannot differentiate LatticeRecord of length $(length(lvr))")
    td = time_domain(lvr)
    LatticeRecord{ET}(lattice(lvr),
        [@. (lvr.records[i + 1] - lvr.records[i])/(td[i + 1] - td[i])
        for i in 1:length(lvr) - 1], (td[2:end] .+ td[1:end-1])./2)
end

"""
    integrate(lattice_record)

Integrates the values stored in the record over time using the trapezoidal rule.
"""
function integrate(lvr::LatticeRecord{ET}) where ET
    td = time_domain(lvr)
    a = copy(first(lvr.records)) * ((td[2] - td[1]) / 2)
    out_records = [zero(first(lvr.records))]
    for i in 2:length(lvr)
        a .+= lvr.records[i] .* (td[i] - td[i-1])
        push!(out_records, @. (a - lvr.records[i] / 2) * (td[i] - td[i-1]))
    end
    LatticeRecord{ET}(lattice(lvr), out_records, td)
end
