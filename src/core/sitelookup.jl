"""
    LookupTable

A helper data structure to quickly find the index of a site in a lattice.

Works well under following assumptions:
- The `sitekey` is some integer property of the sites.
- The sites in the lattice are ordered by `sitekey`.
- The numbering is mostly contiguous, i.e. there are no (or few) gaps in the numbering.
- Secondary keys are also integers, mostly contiguous, ordered and unique for all sites with the same `sitekey`.
Set them to `nothing` to disable usage.
"""
struct LookupTable
    firstkey::Int
    indranges::Vector{UnitRange{Int}}
    secondarykeyranges::Vector{UnitRange{Int}}
end

function sitekey end
secondarykey(::AbstractSite) = nothing
function indrange(lt::LookupTable, site::AbstractSite)
    j = sitekey(site) - lt.firstkey + 1
    1 ≤ j ≤ length(lt.indranges) || return 1:0
    @inbounds irange = lt.indranges[j]
    isempty(lt.secondarykeyranges) && return irange
    @inbounds skrange = lt.secondarykeyranges[j]
    seckey = secondarykey(site)
    start = max(irange.start, irange.stop + seckey - skrange.stop)
    stop = min(irange.stop, seckey - skrange.start + irange.start)
    return start:stop
end

function LookupTable(lat::AbstractLattice)
    isempty(lat) && throw(ArgumentError("Cannot create a lookup table for an empty lattice"))
    firstkey = sitekey(lat[1])
    rowkey = firstkey
    key = rowkey
    rowseckey = secondarykey(lat[1])
    seckey = rowseckey
    skipsec = seckey === nothing

    indranges = Vector{UnitRange{Int}}()
    seckeyranges = Vector{UnitRange{Int}}()
    oldi = 1
    for (i, site) in enumerate(lat)
        key = sitekey(site)
        newseckey = secondarykey(site)
        skipsec = skipsec || newseckey === nothing
        if key != rowkey
            key < rowkey && throw(ArgumentError("The lattice is not ordered by `sitekey`"))
            push!(indranges, oldi:i - 1)
            !skipsec && push!(seckeyranges, rowseckey:seckey)
            for _ in rowkey + 1:key - 1
                push!(indranges, 1:0)
                !skipsec && push!(seckeyranges, 1:0)
            end
            rowkey = key
            rowseckey = newseckey
            oldi = i
        end
        !skipsec && newseckey < rowseckey &&
            throw(ArgumentError("The lattice is not ordered by `secondarykey`"))
        seckey = newseckey
    end
    push!(indranges, oldi:length(lat))
    if skipsec
        empty!(seckeyranges)
    else
        push!(seckeyranges, rowseckey:seckey)
    end
    return LookupTable(firstkey, indranges, seckeyranges)
end

function Base.show(io::IO, ::MIME"text/plain", lt::LookupTable)
    print(io, "Lookup table: ", last(lt.indranges).stop, " sites, ", length(lt.indranges), " batches",
        isempty(lt.secondarykeyranges) ? "" : ", secondary keys enabled")
end

function site_index(lw::LatticeWithParams, site::AbstractSite)
    if hasparam(lw, :lookup)
        lookup_table = getparam(lw, :lookup)
        range = indrange(lookup_table, site)
        isempty(range) && return nothing
        length(range) == 1 && return first(range)
        return site_index(lw.lat, site, range)
    end
    site_index(lw.lat, site)
end

"""
    addlookuptable(lat)

Adds a lookup table to the lattice `lat` and returns the lattice with the lookup table.

!!! Warning
    Make sure you add the lookup table to the lattice after you stop making changes to it.
    Otherwise the results may be unpredictable.

    This operation is not in-place.
"""
function addlookuptable(lat::AbstractLattice)
    return setparam(lat, :lookup, LookupTable(lat))
end
