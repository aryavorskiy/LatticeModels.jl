using LinearAlgebra, Logging, StaticArrays

"""
    Lattice{N,B}
A finite subset of a `B` bravais lattice.

---
    Lattice(sym, sz, bvs[, mask])
Constructs a finite `Lattice{sym, N, NB}` as a subset of the `bvs` Bravais lattice.
`sz` is a `NTuple{N, Int}` which represents how many times the unit cell of `bvs` was translated by each axis - these sites form a *macrocell*.
`mask`, if defined, is a `Vector{Bool}` storing information about which of the sites from the macrocell
are actually included in the lattice, and which are not.

For example, a 3×3 square lattice with its center site excluded is represented as
`Lattice(:square, (3, 3), Bravais([1 0; 0 1]), Bool[1, 1, 1, 1, 0, 1, 1, 1, 1])`

To define a new type of lattice, create an alias for `Lattice{YourSym, YourN, YourNB}`.
Refer to the docs for detailed explanation.
"""
struct Lattice{N, B<:Bravais{Sym, N} where Sym} <: AbstractVector{LatticeSite{N, B}}
    bravais::B
    pointers::Vector{LatticePointer{N}}
end

Base.:(==)(l1::Lattice, l2::Lattice) = (l1.pointers == l2.pointers) && (l1.bravais == l2.bravais)

lattice(l::Lattice) = l
Base.copymutable(l::Lattice) = Lattice(l.bravais, copy(l.pointers))
Base.copy(l::Lattice) = Base.copymutable(l)
macrocell_size(::Lattice) = error("This function is discontinued")
Base.length(l::Lattice) = length(l.pointers)
Base.size(l::Lattice) = (length(l),)
Base.keys(l::Lattice) = Base.OneTo(length(l))
lattice_type(::Lattice{<:Bravais{Sym}}) where {Sym} = Sym
dims(::Lattice{N,B} where B) where {N} = N
dims(l) = dims(lattice(l))
basis_length(l::Lattice) = length(l.bravais)
bravais(l::Lattice) = l.bravais

default_bonds(::Lattice, ::Val) = ()
default_bonds(l::Lattice) = default_bonds(l, Val(1))
default_bonds(l::Lattice, i::Int) = default_bonds(l, Val(i))

get_site(l::Lattice, lp::LatticePointer) = LatticeSite(lp, l.bravais)
get_site(::Lattice, site::LatticeSite) = site
get_site(::Lattice, ::Nothing) = nothing

Base.in(l::Lattice, lp::LatticePointer) = insorted(lp, l.pointers)
Base.in(l::Lattice, site::LatticeSite) = in(l, site.lp)

function Base.getindex(l::Lattice, i::Int)
    @boundscheck checkbounds(l, i)
    return get_site(l, l.pointers[i])
end
Base.getindex(l::Lattice, ci::CartesianIndex{1}) = getindex(l, only(Tuple(ci)))
function Base.getindex(l::Lattice, is::AbstractVector{Int})
    @boundscheck checkbounds(l, is)
    return Lattice(l.bravais, l.pointers[sort(is)])
end
function Base.deleteat!(l::Lattice, inds)
    @boundscheck checkbounds(l, inds)
    deleteat!(l.pointers, inds)
    l
end
Base.pop!(l::Lattice) = Base.deleteat!(l, lastindex(l))
Base.popfirst!(l::Lattice) = Base.deleteat!(l, firstindex(l))

"""
    site_index(l::Lattice, site::LatticeSite; macrocell=false)

Returns the integer index for given `site` in `lattice`.
Returns `nothing` if the site is not present in the lattice.
"""
Base.@propagate_inbounds function site_index(l::Lattice, lp::LatticePointer)
    i = searchsortedfirst(l.pointers, lp)
    i > length(l) && return nothing
    return l.pointers[i] == lp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::Lattice, site::LatticeSite)
    @boundscheck l.bravais == site.bravais
    site_index(l, site.lp)
end
site_index(::Lattice, ::Nothing) = nothing

function Base.iterate(l::Lattice, state = (1, length(l)))
    i, len = state
    return i > len ? nothing : (l[i], (i+1, len))
end

"""
    radius_vector(l::Lattice, site1::LatticeSite, site2::LatticeSite) -> vector
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macrocell to minimize the distance between them).
"""
function radius_vector(l::Lattice, site1::LatticeSite{N}, site2::LatticeSite{N}) where N
    hsz = SVector{N, Int}(macrocell_size(l) .÷ 2)
    tr_unitcell = (site1.unit_cell - site2.unit_cell + hsz) .% macrocell_size(l) - hsz
    bravais(l).basis[:, site1.index] - bravais(l).basis[:, site2.basis_index] + bravais(l).translation_vectors * tr_unitcell
end

"""
    site_distance(l::Lattice, site1::LatticeSite, site2::LatticeSite[; pbc=false])
Returns the distance between two sites on the `l` lattice.

**Keyword arguments:**
- `pbc`: if `true`, the boundary conditions will be considered periodic and
the distance will be measured on the shortest path.
"""
function site_distance(l::Lattice, site1::LatticeSite, site2::LatticeSite; pbc=false)
    if pbc
        norm(radius_vector(l, site1, site2))
    else
        norm(site1.coords - site2.coords)
    end
end

"""
    site_distance(; pbc)
Generates a function that finds the distance between sites (see `site_distance(::Lattice, ::LatticeSite, ::LatticeSite)`).
This notation can be handy when passing this function as an argument.
"""
site_distance(;pbc) = (l, site1, site2) -> site_distance(l, site1, site2, pbc=pbc)

function collect_coords(l::Lattice)
    d = dims(l)
    pts = zeros(d, length(l))
    for (i, site) in enumerate(l)
        pts[:, i] = site.coords
    end
    pts
end

Base.show(io::IO, l::Lattice) = Base.show_default(io, l)
function Base.show(io::IO, ::MIME"text/plain", l::Lattice{N, <:Bravais{Sym}}) where {N,Sym}
    print(io, length(l), "-site ", N, "-dimensional ", Sym, " lattice")
    if basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

"""
    sublattice(lf::Function, l::Lattice) -> Lattice
Generates a a subset of lattice `l` by applying the `lf` function to its sites.
The `lf` function must return a boolean value.
"""
function sublattice(f::Function, l::Lattice)
    inds = [f(site) for site in l]
    return Lattice(l.bravais, l.pointers[inds])
end

function Lattice(bvs::Bravais{Sym,N,NB}, sz::NTuple{N,Int}) where {Sym,N,NB}
    unit_cell = ones(Int, N)
    ptrs = LatticePointer{N}[]
    for unitcell_j in 1:prod(sz)
        svec = SVector{N}(unit_cell)
        for i in 1:NB
            push!(ptrs, LatticePointer(svec, i))
        end
        unitcell_j == prod(sz) && break
        incr_j = N
        while unit_cell[incr_j] == sz[incr_j]
            incr_j -= 1
            incr_j == 0 && break
        end
        unit_cell[incr_j] += 1
        unit_cell[incr_j + 1:end] .= 1
    end
    Lattice(bvs, ptrs)
end

const AnyDimLattice{Sym,NB} = Lattice{N, Bravais{Sym,N,NB}} where N
function Lattice{N, B}(sz::Vararg{Int,N}) where {N,B<:Bravais{Sym,N} where Sym}
    return Lattice(B(), sz)
end
function AnyDimLattice{Sym,NB}(sz::Vararg{Int,N}) where {Sym,N,NB}
    return Lattice(Bravais{Sym,N,NB}(), sz)
end
(::Type{T})(f::Function, sz::Vararg{Int}) where T<:Lattice =
    sublattice(f, T(sz...))

struct IncompatibleLattices <: Exception
    header::String
    l1::Lattice
    l2::Lattice
    IncompatibleLattices(header, l1, l2) = new(header, lattice(l1), lattice(l2))
end
IncompatibleLattices(l1, l2) = IncompatibleLattices("Matching lattices expected", l1, l2)

Base.showerror(io::IO, ex::IncompatibleLattices) = print(io,
"""$(ex.header).\nGot following:
        #1: $(repr("text/plain", ex.l1))
        #2: $(repr("text/plain", ex.l2))""")


"""
Checks if `l1` and `l2` objects are defined on one lattice. Throws an error if not.
"""
function check_samelattice(l1, l2)
    lattice(l1) != lattice(l2) &&
        throw(IncompatibleLattices(l1, l2))
end

"""
Checks if `l1` and `l2` are defined on one macrocell. Throws an error if not.
"""
function check_samebravais(l1, l2)
    la1 = lattice(l1)
    la2 = lattice(l2)
    (la1.bravais != la2.bravais) &&
        throw(IncompatibleLattices("Same lattice type expected", la1, la2))
end

"""
Checks if `l1` is sublattice of `l2`. Throws an error if not.
"""
function check_issublattice(l1::Lattice, l2::Lattice)
    check_samebravais(l1, l2)
    !issubset(l1.pointers, l2.pointers) &&
        throw(IncompatibleLattices("#1 is expected to be sublattice of #2", l1, l2))
end

function Base.union!(l1::Lattice{B, N}, ls::Lattice{B, N}...) where {B, N}
    Base.union!(l1.pointers, (l.pointers for l in ls)...)
    sort!(l1.pointers)
end

function Base.intersect!(l1::Lattice{B, N}, ls::Lattice{B, N}...) where {B, N}
    Base.intersect!(l1.pointers, (l.pointers for l in ls)...)
    sort!(l1.pointers)
end
Base.intersect(l1::Lattice, ls::Lattice...) = Base.intersect!(copy(l1), ls...)

function Base.setdiff!(l1::Lattice{B, N}, ls::Lattice{B, N}...) where {B, N}
    Base.setdiff!(l1.pointers, (l.pointers for l in ls)...)
    sort!(l1.pointers)
end
Base.setdiff(l1::Lattice, ls::Lattice...) = Base.setdiff!(copy(l1), ls...)

Base.emptymutable(l::Lattice{B, N}, ::Type{LatticeSite{N}}=eltype(l)) where {B, N} =
    Lattice(l.bravais, [])
