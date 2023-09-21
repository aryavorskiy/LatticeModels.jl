abstract type AbstractSite{N} end

dims(::AbstractSite{N}) where N = N
Base.iterate(site::AbstractSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)
@inline function Base.getproperty(site::AbstractSite{N}, sym::Symbol) where N
    if sym === :x && N ≥ 1
        return site.coords[1]
    elseif sym === :y && N ≥ 2
        return site.coords[2]
    elseif sym === :z && N ≥ 3
        return site.coords[3]
    end
    Base.getfield(site, sym)
end

abstract type SiteParameter end
get_param(::AbstractSite, p::SiteParameter) = error("Site does not accept param $p")

struct Coord <: SiteParameter axis::Int end
function get_param(site::AbstractSite, c::Coord)
    @assert 1 ≤ c.axis ≤ dims(site)
    return site.coords[c.axis]
end
SiteParameter(sym::Symbol) =
    sym === :x ? Coord(1) : sym === :y ? Coord(2) : sym === :z ? Coord(3) :
        error("Failed to parse parameter `$sym`. Try using `p\"$sym\"`.")

abstract type AbstractLattice{SiteT} <: AbstractSet{SiteT} end
lattice(l::AbstractLattice) = l
dims(::AbstractLattice{<:AbstractSite{N}}) where {N} = N
dims(l) = dims(lattice(l))
Base.size(l::AbstractLattice) = (length(l),)
site_index(::AbstractLattice, ::Nothing) = nothing

# Set functions
Base.copy(l::AbstractLattice) = Base.copymutable(l)
Base.in(site::SiteT, l::AbstractLattice{SiteT}) where {SiteT} =
    site_index(l, site) !== nothing
Base.hasfastin(::Type{<:AbstractLattice}) = true
function Base.delete!(l::AbstractLattice{ST}, site::ST) where ST<:AbstractSite
    i = site_index(l, site)
    i !== nothing && Base.deleteat!(l, i)
    return l
end

# Indexing
Base.firstindex(::AbstractLattice) = 1
Base.lastindex(l::AbstractLattice) = length(l)
Base.eachindex(l::AbstractLattice) = firstindex(l):lastindex(l)
Base.checkbounds(::Type{Bool}, l::AbstractLattice, is) = all(1 .≤ is .≤ length(l))
Base.checkbounds(l::AbstractLattice, is) =
    !Base.checkbounds(Bool, l, is) && throw(BoundsError(l, is))
Base.pop!(l::AbstractLattice) = Base.deleteat!(l, lastindex(l))
Base.popfirst!(l::AbstractLattice) = Base.deleteat!(l, firstindex(l))

function pairs_to_inds(l::AbstractLattice, pairs...)
    inds = Int[]
    for (i, site) in enumerate(l)
        if all(get_param(site, param) in val for (param, val) in pairs)
            push!(inds, i)
        end
    end
    return inds
end
function Base.getindex(l::AbstractLattice, pairs::Pair{<:SiteParameter}...)
    inds = pairs_to_inds(l, pairs...)
    return length(inds) == 1 ? l[only(inds)] : l[inds]
end
function Base.getindex(l::AbstractLattice; kw...)
    return l[(SiteParameter(crd) => val for (crd, val) in kw)...]
end

function collect_coords(l::AbstractLattice)
    d = dims(l)
    pts = zeros(d, length(l))
    for (i, site) in enumerate(l)
        pts[:, i] = site.coords
    end
    pts
end

"""
    sublattice(lf::Function, l::Lattice) -> Lattice
Generates a a subset of lattice `l` by applying the `lf` function to its sites.
The `lf` function must return a boolean value.
"""
function sublattice(f::Function, l::AbstractLattice)
    inds = [i for (i, site) in enumerate(l) if f(site)]
    return l[inds]
end

struct IncompatibleLattices <: Exception
    header::String
    l1::AbstractLattice
    l2::AbstractLattice
    IncompatibleLattices(header, l1, l2) = new(header, lattice(l1), lattice(l2))
end
IncompatibleLattices(l1, l2) = IncompatibleLattices("Incompatible lattices", l1, l2)

Base.showerror(io::IO, ex::IncompatibleLattices) = print(io,
"""$(ex.header).\nGot following:
        #1: $(repr("text/plain", ex.l1))
        #2: $(repr("text/plain", ex.l2))""")

"""
Checks if `l1` and `l2` objects are defined on one lattice. Throws an error if not.
"""
function check_samesites(l1, l2)
    sites(lattice(l1)) != sites(lattice(l2)) &&
        throw(IncompatibleLattices("Matching lattices expected", l1, l2))
end

"""
Checks if `l1` is sublattice of `l2`. Throws an error if not.
"""
function check_issublattice(l1::AbstractLattice, l2::AbstractLattice)
    !issubset(l1, l2) &&
        throw(IncompatibleLattices("#1 is expected to be sublattice of #2", l1, l2))
end

"""
    site_distance(l::Lattice, site1::LatticeSite, site2::LatticeSite[; pbc=false])
Returns the distance between two sites on the `l` lattice.

**Keyword arguments:**
- `pbc`: if `true`, the boundary conditions will be considered periodic and
the distance will be measured on the shortest path.
"""
function site_distance(site1::AbstractSite, site2::AbstractSite)
    norm(site1.coords - site2.coords)
end

"""
    site_distance(; pbc)
Generates a function that finds the distance between sites (see `site_distance(::Lattice, ::LatticeSite, ::LatticeSite)`).
This notation can be handy when passing this function as an argument.
"""
site_distance() = (site1, site2) -> site_distance(site1, site2)

shift_site(::AbstractLattice, site::AbstractSite) = (1, site)
const SingleBond{LT<:AbstractSite} = Pair{LT, LT}
