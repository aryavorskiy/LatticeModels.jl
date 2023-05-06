using RecipesBase, LinearAlgebra, Logging, StaticArrays
import Base: length, size, copy, iterate, getindex, eltype, show, ==, isless,
    pop!, popat!, popfirst!, splice!, lastindex, getproperty

"""
    Bravais{N, NB}
`N`-dimensional infinite Bravais lattice with `NB` sites in basis.

---
    Bravais(translation_vectors[, basis])
Constructs a Bravais lattice with given translation vectors and locations of basis sites
relative to some unit cell.
The `basis` argument can be omitted, in which case the lattice basis will consist of
one site located in the bottom-left corner of the unit cell.

`translation_vectors` argument must be an `AbstractMatrix{<:Real}` of size `N×N`,
while `basis` must also be an  abstract matrix of size `N×NB`.
"""
struct Bravais{N,NB}
    translation_vectors::SMatrix{N,N,Float64}
    basis::SMatrix{N,NB,Float64}
    function Bravais(translation_vectors::AbstractMatrix{<:Real}, basis::AbstractMatrix{<:Real},
        origin::AbstractVector{<:Real}=zeros(size(basis)[1]))
        (size(translation_vectors)[1] != size(basis)[1]) &&
            error("inconsistent dimension count (got $(size(translation_vectors)[1]), $(size(basis)[1]))")
        N, NB = size(basis)
        new{N,NB}(translation_vectors, basis .- origin)
    end
end
Bravais(translation_vectors::AbstractMatrix{<:Real}) =
    Bravais(translation_vectors, zeros((size(translation_vectors)[1], 1)))

dims(@nospecialize _::Bravais{N}) where {N} = N
length(::Bravais{N,NB}) where {N,NB} = NB

"""
    Lattice{LatticeSym, N, NB}
A finite subset of a `Brvais{N, NB}`. `LatticeSym` is a `Symbol` which represents
the type of the lattice (e. g. `:square`, `:honeycomb`).
This makes `Lattice` object behavior known at compile-time,
which allows to introduce various optimizations or to define specific plot recipes.

---
    Lattice(sym, sz, bvs[, mask])
Constructs a finite `Lattice{sym, N, NB}` as a subset of the `bvs` Bravais lattice.
`sz` is a `NTuple{N, Int}` which represents how many times the unit cell of `bvs` was translated by each axis - these sites form a *macro cell*.
`mask`, if defined, is a `Vector{Bool}` storing information about which of the sites from the macro cell
are actually included in the lattice, and which are not.

For example, a 3×3 square lattice with its center site excluded is represented as
`Lattice(:square, (3, 3), Bravais([1 0; 0 1]), Bool[1, 1, 1, 1, 0, 1, 1, 1, 1])`

To define a new type of lattice, create an alias for `Lattice{YourSym, YourN, YourNB}`.
Refer to the docs for detailed explanation.
"""
struct Lattice{LatticeSym,N,NB}
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

==(l1::T, l2::T) where {T<:Lattice} =
    (size(l1) == size(l2)) && (bravais(l1) == bravais(l2))

copy(l::Lattice{LatticeSym}) where {LatticeSym} =
    Lattice(LatticeSym, size(l), bravais(l), copy(l.mask))
size(l::Lattice) = l.lattice_size
length(l::Lattice) = count(l.mask)
lattice_type(::Lattice{LatticeSym}) where {LatticeSym} = LatticeSym
dims(::Lattice{LatticeSym,N} where {LatticeSym}) where {N} = N
dims(l) = dims(lattice(l))
basis_length(::Lattice{LatticeSym,N,NB} where {LatticeSym,N}) where {NB} = NB
bravais(l::Lattice) = l.bravais

site_coords(l::Lattice, basis_index::Int, unit_cell) =
    bravais(l).basis[:, basis_index] + bravais(l).translation_vectors * unit_cell
site_coords(l::Lattice{Sym,N,1} where {Sym,N}, ::Int, unit_cell) =
    bravais(l).translation_vectors * unit_cell
site_coords(::Lattice{:square,N,1} where {N}, ::Int, unit_cell) = Float64.(unit_cell)

"""
    LatticeSite{N}
A site of a `Lattice{LatticeSym, N, NB}` lattice.

Fields:
- `unit_cell`: a set of translations along all axes representing the unit cell the site is located in.
- `basis_index`: the number of site in the lattice basis.

This type is used to iterate over all sites of a `Lattice{LatticeSym, N, NB}`.
The exact location of a `LatticeSite` can be found using the `site.coords` function.
"""
struct LatticeSite{N}
    unit_cell::SVector{N,Int}
    basis_index::Int
    coords::SVector{N,Float64}
end
function getproperty(site::LatticeSite{N}, sym::Symbol) where N
    if N ≥ 1 && sym === :x
        site.coords[1]
    elseif N ≥ 2 && sym === :y
        site.coords[2]
    elseif N ≥ 3 && sym === :z
        site.coords[3]
    else
        getfield(site, sym)
    end
end

function LatticeSite(unit_cell, basis_index, l::Lattice)
    LatticeSite(unit_cell, basis_index, site_coords(l, basis_index, unit_cell))
end

iterate(site::LatticeSite{N}, i=1) where N = i > N ? nothing : (site.coords[i], i + 1)

_cind(site::LatticeSite) = CartesianIndex(site.basis_index, site.unit_cell...)

_cinds(l::Lattice{LatticeSym,N,NB} where {LatticeSym}) where {N,NB} =
    CartesianIndex{N + 1}(1):CartesianIndex(NB, size(l)...)

site_coords(::Lattice, site::LatticeSite) =  error("`site_coords(l, site)` is no longer available. Use `site.coords` instead")

function getindex(l::Lattice{Sym, N} where Sym, i::Int) where N
    counter = 0
    cinds = _cinds(l)
    i ≤ 0 && throw(BoundsError(l, i))
    for j in 1:length(l.mask)
        counter += l.mask[j]
        if counter == i
            basis_index, unit_cell... = Tuple(cinds[j])
            return LatticeSite(SVector(unit_cell), basis_index, l)
        end
    end
    throw(BoundsError(l, i))
end
getindex(l::Lattice, ci::CartesianIndex{1}) = getindex(l, only(Tuple(ci)))
function getindex(l::Lattice{Sym}, is::AbstractVector{Int}) where Sym
    new_mask = zero(l.mask)
    try
        (@view new_mask[l.mask])[is] .= true
    catch BoundsError
        throw(BoundsError(l, is))
    end
    Lattice(Sym, size(l), bravais(l), new_mask)
end

"""
    site_index(l::Lattice, site::LatticeSite; macrocell=false)

Returns the integer index for given `site` in `lattice`.
Returns `nothing` if the site is not present in the lattice.
"""
function site_index(l::Lattice, site::LatticeSite)
    i = LinearIndices((1:basis_length(l), (1:s for s in size(l))...))[_cind(site)]
    (i === nothing || !l.mask[i]) && return nothing
    count(@view l.mask[1:i])
end
site_index(::LatticeSite, ::Lattice) = error("site_index(::LatticeSite, ::Lattice) discontinued." *
"This error message will be removed in future versions. Use site_index(::Lattice, ::LatticeSite) instead.")

function splice!(l::Lattice, is)
    view(l.mask, l.mask)[collect(is)] .= false
    l
end
popat!(l::Lattice, i::Int) = splice!(l, i)
pop!(l::Lattice) = popat!(l, length(l))
popfirst!(l::Lattice) = popat!(l, 1)
lastindex(l::Lattice) = length(l)

Base.eltype(::Lattice{LatticeSym,N}) where {LatticeSym,N} = LatticeSite{N}

==(@nospecialize(site1::LatticeSite), @nospecialize(site2::LatticeSite)) =
    site1.basis_index == site2.basis_index && site1.unit_cell == site2.unit_cell

isless(@nospecialize(site1::LatticeSite), @nospecialize(site2::LatticeSite)) =
    isless(_cind(site1), _cind(site2))

function iterate(l::Lattice{Sym,N} where Sym) where N
    cinds = _cinds(l)
    index = findfirst(l.mask)
    index === nothing && return nothing
    cind = cinds[index]
    basis_index, unit_cell... = Tuple(cind)
    LatticeSite(SVector(unit_cell), basis_index, l), (cinds, index)
end

function iterate(l::Lattice{Sym,N} where Sym, state) where N
    cinds, index = state
    index = findnext(l.mask, index + 1)
    index === nothing && return nothing
    cind = cinds[index]
    basis_index, unit_cell... = Tuple(cind)
    LatticeSite(SVector(unit_cell), basis_index, l), (cinds, index)
end

"""
    radius_vector(l::Lattice, site1::LatticeSite, site2::LatticeSite) -> vector
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macro cell to minimize the distance between them).
"""
function radius_vector(l::Lattice, site1::LatticeSite{N}, site2::LatticeSite{N}) where N
    hsz = SVector{N, Int}(size(l) .÷ 2)
    tr_unitcell = (site1.unit_cell - site2.unit_cell + hsz) .% size(l) - hsz
    bravais(l).basis[:, site1.basis_index] - bravais(l).basis[:, site2.basis_index] + bravais(l).translation_vectors * tr_unitcell
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

@recipe function f(l::Lattice; pretty=true, high_contrast=false)
    if high_contrast
        pretty = false
        markersize := 4
        markercolor := :black
        markerstrokealpha := 1
        markerstrokestyle := :solid
        markerstrokewidth := 2
        markerstrokecolor := :white
    end
    if pretty
        l_outp = copy(l)
        fill!(l_outp.mask, true)
        annotations = repeat(Any[""], length(l_outp))
        annotations[l.mask] .= ((i, :left, :top, :grey, 8) for i in 1:length(l))
        series_annotations := annotations
        seriesalpha := l.mask .* 0.9 .+ 0.1
        label --> ""
        l_outp, nothing
    else
        l, nothing
    end
end

function collect_coords(l::Lattice)
    d = dims(l)
    pts = zeros(d, length(l))
    for (i, site) in enumerate(l)
        pts[:, i] = site.coords
    end
    pts
end

@recipe function f(l::Lattice, v)
    aspect_ratio := :equal
    marker_z := v
    pts = collect_coords(l)
    if dims(l) == 3
        X, Y, Z = eachrow(pts)
        Xr, Yr, Zr = eachrow(round.(pts, digits=3))
        seriestype := :scatter3d
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ", ", Zr, ")")
        end
        X, Y, Z
    else
        if dims(l) == 1
            X = vec(pts)
            Y = zero(X)
        else
            X, Y = eachrow(pts[1:2, :])
        end
        seriestype --> :scatter
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            Xr, Yr = eachrow(round.(pts, digits=3))
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ")")
        end
        if plotattributes[:seriestype] == :scatter
            X, Y
        elseif plotattributes[:seriestype] == :surface
            X, Y, v
        else
            throw(ArgumentError("unsupported series type $(plotattributes[:seriestype])"))
        end
    end
end

function show(io::IO, ::MIME"text/plain", l::Lattice{LatticeSym,N}) where {N,LatticeSym}
    print(io, "$(length(l))-site ", LatticeSym)
    if N != 1
        print(io, " lattice on ", join(size(l), "×"), " macro cell")
    elseif !all(l.mask)
        print(io, " chain on ", size(l)[1], " unit cells")
    else
        print(io, " chain")
    end
    if N > 1 && basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

"""
    sublattice(lf::Function, l::Lattice) -> Lattice
Generates a a subset of lattice `l` by applying the `lf` function to its sites.
The `lf` function must return a boolean value.
"""
function sublattice(f::Function, l::Lattice{LatticeSym}) where {LatticeSym}
    new_mask = zero(l.mask)
    new_mask[l.mask] = [f(site) for site in l]
    Lattice(LatticeSym, size(l), bravais(l), new_mask)
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

function check_lattice_match(l1, l2)
    lattice(l1) != lattice(l2) &&
        throw(ArgumentError("""lattice mismatch:
        $(repr("text/plain", lattice(l1)))
        $(repr("text/plain", lattice(l2)))"""))
end

function check_macrocell_match(l1, l2)
    la1 = lattice(l1)
    la2 = lattice(l2)
    (size(la1) != size(la2) || bravais(la1) != bravais(la2)) &&
        throw(ArgumentError("""macrocell mismatch:
        $(size(la1))-size with $(bravais(la1))
        $(size(la2))-size with $(bravais(la2))"""))
end

function check_is_sublattice(l1::Lattice, l2::Lattice)
    check_macrocell_match(l1, l2)
    any(.!l1.mask .& l2.mask) &&
        error("macrocells match but sublattice check failed")
end
