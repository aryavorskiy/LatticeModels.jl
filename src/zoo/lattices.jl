"""
    Chain(sz)

Construct a 1D chain lattice of size `sz`.
"""
Chain(f::Function, axes; kw...) = span_unitcells(f, UnitCell(SMatrix{1,1}(1.0)), axes; kw...)
Chain(axes; kw...) = Chain(alwaystrue, axes; kw...)

abstract type BravaisLatticeType end
abstract type BravaisLatticeTypeVarN{N} <: BravaisLatticeType end
function construct_unitcell(T::Type{<:BravaisLatticeType}, NU)
    uc = construct_unitcell(T)
    if ldims(uc) != NU
        throw(ArgumentError("Invalid dimensionality for $T; expected $(ldims(uc))"))
    end
    return uc
end
_try_apply_typevar(T::Type{<:BravaisLatticeTypeVarN}, NU) = T{NU}
_try_apply_typevar(T::Type{<:BravaisLatticeTypeVarN{NU}}, _NU) where NU = T
construct_unitcell(T::Type{<:BravaisLatticeTypeVarN}, NU) =
    construct_unitcell(_try_apply_typevar(T, NU))

"""
    @bravaisdef MyBravaisLattice UnitCell(...)
    @bravaisdef MyBravaisLattice N -> UnitCell(...)

Define a new Bravais lattice type `MyBravaisLattice` with a unit cell constructor `UnitCell(expr)`.
If the notation is `N -> UnitCell(expr)`, the unit cell constructor will be dependent on the dimensionality `N`.
otherwise, the dimensionality will be inferred from the unit cell.
`N` is the dimensionality of the lattice.

## Examples

```jldoctest
julia> using LatticeModels

julia> @bravaisdef MyBravaisLattice UnitCell([1 0; 0 1]);   # 2D square lattice

julia> MyBravaisLattice(3, 3)
9-site 2-dim Bravais lattice in 2D space
Unit cell:
  Basis site coordinates:
    ┌      ┐
    │ 0.000│
    │ 0.000│
    └      ┘
  Translation vectors:
    ┌      ┐ ┌      ┐
    │ 1.000│ │ 0.000│
    │ 0.000│ │ 1.000│
    └      ┘ └      ┘
Lattice type: MyBravaisLattice
Default translations:
  :axis1 → Bravais[3, 0]
  :axis2 → Bravais[0, 3]
Nearest neighbor hoppings:
  1.00000 =>
    Bravais[1, 0]
    Bravais[0, 1]
  1.41421 =>
    Bravais[1, -1]
    Bravais[1, 1]
  2.00000 =>
    Bravais[2, 0]
    Bravais[0, 2]
Boundary conditions: none
```
"""
macro bravaisdef(type, expr)
    is_ndep = Meta.isexpr(expr, :->)
    if is_ndep
        @assert expr.args[1] isa Symbol "Invalid unitcell constructor; one parameter N expected"
        N_sym = expr.args[1]
        type_expr = :($type{$N_sym})
        unitcell_construct = expr.args[2]
        return quote
            Core.@__doc__ struct $type_expr <: LatticeModels.BravaisLatticeTypeVarN{$N_sym} end
            LatticeModels.construct_unitcell(::Type{$type_expr}) where $N_sym =
                $unitcell_construct
            LatticeModels.construct_unitcell(::Type{$type}) =
                throw(ArgumentError("Unknown dimensionality for " * $(string(type)) *
                    "; Try using " * $(string(type)) * "{N}"))
            LatticeModels.latticename(::Type{$type_expr}) where $N_sym = $(string(type))
        end |> esc
    else
        return quote
            Core.@__doc__ struct $type <: LatticeModels.BravaisLatticeType end
            LatticeModels.construct_unitcell(::Type{$type}) = $expr
            LatticeModels.latticename(::Type{$type}) = $(string(type))
        end |> esc
    end
end
UnitCell(::Type{T}; kw...) where T<:BravaisLatticeType = transform_unitcell(construct_unitcell(T); kw...)
function (::Type{T})(f::Function, sz::Vararg{LatticeModels.RangeT, NU}; kw...) where {T<:BravaisLatticeType,NU}
    l = LatticeModels.span_unitcells(f, construct_unitcell(T, NU), sz...; kw...)
    if T <: BravaisLatticeTypeVarN
        return settype(l, _try_apply_typevar(T, NU))
    else
        return settype(l, T)
    end
end
function (::Type{T})(args...; kw...) where T<:BravaisLatticeType
    T(alwaystrue, args...; kw...)
end

Base.show(io::IO, ::MIME"text/plain", T::Type{<:BravaisLatticeType}) = print(io, "Lattice type: ", string(T))
settype(l::AbstractLattice, T::Type{<:BravaisLatticeType}) = pushmeta(l, :latticetype, T)
gettype(l::AbstractLattice) = getmeta(l, :latticetype, nothing)
function checktype(l::AbstractLattice, T::Type{<:BravaisLatticeType})
    AT = gettype(l)
    if AT === nothing
        @warn "Lattice type not defined; expected $T"
    elseif !(AT <: T)
        throw(ArgumentError("Invalid lattice type $AT; expected $T"))
    end
end
checktype(any, T) = checktype(lattice(any), T)

function _ngonshape(r::SVector{2}, n, skip=true, expand=true)
    expand && (r = r * 1.04)
    is = (skip / 2:n + skip / 2) .* 2pi/n
    skip && (r /= cos(pi / n))
    return r[1] .* cos.(is) .+ r[2] .* sin.(is), r[1] .* sin.(is) .- r[2] .* cos.(is)
end

"""
    SquareLattice{N}
Represents a square lattice in `N` dimensions.

---
    SquareLattice(sz...)

Construct a square lattice of size `sz`.
"""
@bravaisdef SquareLattice N -> UnitCell(SMatrix{N,N}(I))
getshape(l::BravaisLattice, ::Type{SquareLattice{2}}) =
    _ngonshape(unitvector(l, 1) / 2, 4)

"""
    TriangularLattice
Represents a triangular lattice.
Lattice vectors: `[1, 0]` and `[0.5, √3/2]`.

---
    TriangularLattice(a, b)

Construct a triangular lattice of a×b spanned unit cells.
"""
@bravaisdef TriangularLattice UnitCell([1 0.5; 0 √3/2])
getshape(l::BravaisLattice, ::Type{TriangularLattice}) =
    _ngonshape(unitvector(l, 1) / 2, 6)

"""
    HoneycombLattice
Represents a honeycomb lattice.

Lattice vectors: `[1, 0]` and `[0.5, √3/2]`,
two sites at `[0, 0]` and `[0.5, √3/6]` in each unit cell.

---
    HoneycombLattice(a, b)

Construct a honeycomb lattice of a×b spanned unit cells.
"""
@bravaisdef HoneycombLattice UnitCell([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
getshape(l::BravaisLattice, ::Type{HoneycombLattice}) =
    _ngonshape(unitvector(l, 1) / 3, 6, false)

"""
    GrapheneRibbon(len, wid[, center; kw...])

Construct a graphene ribbon sample with zigzag edges.
To get armchair edges, simply rotate the lattice by 90 degrees.

## Arguments
- `len`: the length of the ribbon.
- `wid`: the width of the ribbon.
- `center`: the unit cell coordinates of the bottom-left corner of the ribbon. Default is `(0, 0)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function GrapheneRibbon(len, wid, center=(0, 0); kw...)
    j1ind = center[1] - wid:len + center[1] - 1
    j2ind = center[2]:wid + center[2] - 1
    default_translations = wid % 2 == 0 ?
        (:horizontal => Bravais[len, 0], :vertical => Bravais[-wid ÷ 2, wid]) :
        (:horizontal => Bravais[len, 0])
    return HoneycombLattice(j1ind, j2ind; unitvectortrs=false,
            default_translations=default_translations, kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2 / 2
        if -0.5 ≤ j_prj ≤ len - 0.5
            return !(site.basindex == 2 && j_prj == len - 0.5 ||
                site.basindex == 1 && j_prj == -0.5)
        else
            return false
        end
    end
end

"""
    KagomeLattice
Represents a kagome lattice.

Lattice vectors: `[1, 0]` and `[0.5, √3/2]`,
three sites at `[0, 0]`, `[0.5, 0]` and `[0.25, √3/4]` in each unit cell.

---
    KagomeLattice(a, b)

Construct a kagome lattice of a×b spanned unit cells.
"""
@bravaisdef KagomeLattice UnitCell([1 0.5; 0 √3/2], [0 0.5 0.25; 0 0 √3/4])
getshape(l::BravaisLattice, ::Type{KagomeLattice}) =
    _ngonshape(unitvector(l, 1) / 4, 6)
