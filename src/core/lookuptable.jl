"""
    LookupTable

A helper data structure to quickly find the index of a site in a lattice.

Relies on `sitekey(site)` and `secondarykey(function)` functions to determine the index of a site in the lattice.

Works well under following assumptions:
- The `sitekey` is some integer property of the sites.
- The sites in the lattice are ordered by `sitekey`.
- The numbering is mostly contiguous, i.e. there are no (or few) gaps in the numbering.
- The `secondarykey` is also integer, mostly contiguous, ordered and unique for all sites with the same `sitekey`.
Set them to `nothing` to disable usage.
"""
struct LookupTable
    firstkey::Int
    strides::Vector{UnitRange{Int}}
    secondarykeyranges::Vector{UnitRange{Int}}
end

sitekey(site::AbstractSite) =
    throw(ArgumentError("Cannot build a lookup table for site type `$(typeof(site))`"))
secondarykey(::AbstractSite) = nothing
function indrange(lt::LookupTable, site::AbstractSite)
    j = sitekey(site) - lt.firstkey + 1
    1 ≤ j ≤ length(lt.strides) || return 1:0
    @inbounds range = lt.strides[j]
    isempty(lt.secondarykeyranges) && return range
    @inbounds skrange = lt.secondarykeyranges[j]
    seckey = secondarykey(site)
    start = max(range.start, range.stop + seckey - skrange.stop)
    stop = min(range.stop, seckey - skrange.start + range.start)
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
    print(io, "Lookup table: ", last(lt.strides).stop, " sites, ", length(lt.strides), " strides",
        isempty(lt.secondarykeyranges) ? "" : ", secondary keys enabled")
end

function site_index(lw::LatticeWithMetadata, site::AbstractSite, range)
    if hasmeta(lw, :lookup)
        lookup_table = getmeta(lw, :lookup)
        lrange = indrange(lookup_table, site)
        tot_range = intersect(range, lrange)
        isempty(tot_range) && return nothing
        length(tot_range) == 1 && return first(tot_range)
        return site_index(lw.lat, site, tot_range)
    end
    site_index(lw.lat, site)
end

"""
    addlookuptable(lat)

Adds a lookup table to the lattice `lat` and returns the lattice with the lookup table.

!!! warning
    Make sure you add the lookup table to the lattice after you stop making changes to it.
    Otherwise the results may be unpredictable.

    This operation is not in-place.
"""
function addlookuptable(lat::AbstractLattice)
    return setmeta(lat, :lookup, LookupTable(lat))
end
