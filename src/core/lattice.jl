using LinearAlgebra, Logging, StaticArrays

"""
    Lattice{LatticeSym, N, NB}
A finite subset of a `Brvais{N, NB}`. `LatticeSym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Lattice` object behavior known at compile-time,
which allows to introduce various optimizations or to define specific plot recipes.

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
struct Lattice{LatticeSym,N,NB} <: AbstractVector{LatticeSite{N}}
    lattice_size::NTuple{N,Int}
    bravais::Bravais{N,NB}
    mask::Vector{Bool}
    function Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N,NB}, mask; kw...) where {N,NB}
        length(mask) != prod(sz) * length(bvs) &&
            error("inconsistent mask length")
        new{sym,N,NB}(sz, bvs, mask)
    end
end
Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}) where {N} =
    Lattice(sym::Symbol, sz::NTuple{N,Int}, bvs::Bravais{N}, fill(true, prod(sz) * length(bvs)))

Base.:(==)(l1::T, l2::T) where {T<:Lattice} =
    (macrocell_size(l1) == macrocell_size(l2)) && (bravais(l1) == bravais(l2))

lattice(l::Lattice) = l
Base.copymutable(l::Lattice{LatticeSym}) where {LatticeSym} =
    Lattice(LatticeSym, macrocell_size(l), bravais(l), copy(l.mask))
Base.copy(l::Lattice) = Base.copymutable(l)
macrocell_size(l::Lattice) = l.lattice_size
Base.length(l::Lattice) = count(l.mask)
Base.size(l::Lattice) = (count(l.mask),)
Base.keys(l::Lattice) = Base.OneTo(length(l))
lattice_type(::Lattice{LatticeSym}) where {LatticeSym} = LatticeSym
dims(::Lattice{LatticeSym,N} where {LatticeSym}) where {N} = N
dims(l) = dims(lattice(l))
basis_length(::Lattice{LatticeSym,N,NB} where {LatticeSym,N}) where {NB} = NB
bravais(l::Lattice) = l.bravais

site_coords(l::Lattice, lp::LatticePointer) =
    bravais(l).basis[:, lp.basis_index] + bravais(l).translation_vectors * lp.unit_cell
site_coords(l::Lattice{Sym,N,1} where {Sym,N}, lp::LatticePointer) =
    bravais(l).translation_vectors * lp.unit_cell
site_coords(::Lattice{:square,N,1} where {N}, lp::LatticePointer) = Float64.(lp.unit_cell)

default_bonds(::Lattice, ::Val) = ()
default_bonds(l::Lattice) = default_bonds(l, Val(1))

get_site(l::Lattice, lp::LatticePointer) = LatticeSite(lp, site_coords(l, lp))
get_site(::Lattice, site::LatticeSite) = site
get_site(::Lattice, ::Nothing) = nothing

cartesian_indices(l::Lattice{LatticeSym,N,NB} where {LatticeSym}) where {N,NB} =
    CartesianIndex{N + 1}(1):CartesianIndex(NB, macrocell_size(l)...)
linear_indices(l::Lattice{LatticeSym,N,NB} where {LatticeSym,N}) where {NB} =
    LinearIndices((NB, macrocell_size(l)...))

site_coords(::Lattice, site::LatticeSite) =  error("`site_coords(l, site)` is no longer available. Use `site.coords` instead")

function Base.getindex(l::Lattice{Sym, N} where Sym, i::Int) where N
    counter = 0
    cinds = cartesian_indices(l)
    i ≤ 0 && throw(BoundsError(l, i))
    for j in 1:length(l.mask)
        counter += l.mask[j]
        if counter == i
            return get_site(l, LatticePointer(cinds[j]))
        end
    end
    throw(BoundsError(l, i))
end
Base.getindex(l::Lattice, ci::CartesianIndex{1}) = getindex(l, only(Tuple(ci)))
function Base.getindex(l::Lattice{Sym}, is::AbstractVector{Int}) where Sym
    new_mask = zero(l.mask)
    @boundscheck checkbounds(l, is)
    (@view new_mask[l.mask])[is] .= true
    Lattice(Sym, macrocell_size(l), bravais(l), new_mask)
end
function Base.deleteat!(l::Lattice, inds)
    @boundscheck checkbounds(l, inds)
    view(l.mask, l.mask)[collect(inds)] .= false
    l
end
Base.pop!(l::Lattice) = Base.deleteat!(l, lastindex(l))
Base.popfirst!(l::Lattice) = Base.deleteat!(l, firstindex(l))

"""
    site_index(l::Lattice, site::LatticeSite; macrocell=false)

Returns the integer index for given `site` in `lattice`.
Returns `nothing` if the site is not present in the lattice.
"""
function site_index(l::Lattice, site::Union{LatticeSite, LatticePointer})
    linds = linear_indices(l)
    cind = cartesian_index(site)
    i = get(linds, cind, nothing)
    (i === nothing || !l.mask[i]) && return nothing
    count(@view l.mask[1:i])
end
site_index(::Lattice, ::Nothing) = nothing

function Base.iterate(l::Lattice{Sym,N} where Sym) where N
    cinds = cartesian_indices(l)
    index = findfirst(l.mask)
    index === nothing && return nothing
    get_site(l, LatticePointer(cinds[index])), (cinds, index)
end

function Base.iterate(l::Lattice{Sym,N} where Sym, state) where N
    cinds, index = state
    index = findnext(l.mask, index + 1)
    index === nothing && return nothing
    get_site(l, LatticePointer(cinds[index])), (cinds, index)
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
function Base.show(io::IO, ::MIME"text/plain", l::Lattice{LatticeSym,N}) where {N,LatticeSym}
    print(io, "$(length(l))-site ", LatticeSym)
    if N != 1
        print(io, " lattice (", join(macrocell_size(l), "×"), " macrocell")
    elseif !all(l.mask)
        print(io, " chain (", macrocell_size(l)[1], " unit cells")
    else
        print(io, " chain")
    end
    if N > 1 && basis_length(l) > 1
        print(io, ", ", basis_length(l), "-site basis")
    end
    print(io, ")")
end

"""
    sublattice(lf::Function, l::Lattice) -> Lattice
Generates a a subset of lattice `l` by applying the `lf` function to its sites.
The `lf` function must return a boolean value.
"""
function sublattice(f::Function, l::Lattice{LatticeSym}) where {LatticeSym}
    new_mask = zero(l.mask)
    new_mask[l.mask] = [f(site) for site in l]
    Lattice(LatticeSym, macrocell_size(l), bravais(l), new_mask)
end

function Lattice{T,N,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

const AnyDimLattice{T,NB} = Lattice{T,N,NB} where {N}

AnyDimLattice{T,NB}(sz::Vararg{Int,N}; kw...) where {T,N,NB} =
    Lattice{T,N,NB}(sz...; kw...)

function AnyDimLattice{T,NB}(f::Function, sz::Vararg{Int,N}; kw...) where {T,N,NB}
    l = Lattice{T,N,NB}(sz...; kw...)
    sublattice(f, l)
end

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
function check_samemacrocell(l1, l2)
    la1 = lattice(l1)
    la2 = lattice(l2)
    (macrocell_size(la1) != macrocell_size(la2) || bravais(la1) != bravais(la2)) &&
        throw(IncompatibleLattices("Lattices on matching macrocell expected", la1, la2))
end

"""
Checks if `l1` is sublattice of `l2`. Throws an error if not.
"""
function check_issublattice(l1::Lattice, l2::Lattice)
    check_samemacrocell(l1, l2)
    any(l1.mask .& .!l2.mask) &&
        throw(IncompatibleLattices("#1 is expected to be sublattice of #2", l1, l2))
end

function Base.union!(l1::Lattice{Sym}, l2::Lattice{Sym}) where Sym
    check_samemacrocell(l1, l2)
    @. l1.mask = l1.mask | l2.mask
    l1
end

function Base.intersect!(l1::Lattice{Sym}, l2::Lattice{Sym}) where Sym
    check_samemacrocell(l1, l2)
    @. l1.mask = l1.mask & l2.mask
    l1
end
Base.intersect(l1::Lattice{Sym}, l2::Lattice{Sym}) where Sym =
    Base.intersect!(copy(l1), l2)

function Base.setdiff!(l1::Lattice{Sym}, l2::Lattice{Sym}) where Sym
    check_samemacrocell(l1, l2)
    @. l1.mask = l1.mask & !l2.mask
    l1
end

Base.emptymutable(l::Lattice{Sym, N}, ::Type{LatticeSite{N}}=eltype(l)) where {Sym, N} =
    Lattice(Sym, macrocell_size(l), bravais(l), zero(l.mask))
