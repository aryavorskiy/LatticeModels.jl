function construct_unitcell end

function _add_parameter_to_call!(expr, sym)
    if Meta.isexpr(expr, (:block, :if, :elseif, :return))
        foreach(e -> _add_parameter_to_call!(e, sym), expr.args)
    elseif Meta.isexpr(expr, :call)
        if expr.args[1] == :UnitCell
            expr.args[1] = :(LatticeModels.UnitCell{$(QuoteNode(sym))})
        end
    end
end

"""
    @bravaisdef MyBravaisLattice UnitCell(...)
    @bravaisdef MyBravaisLattice N -> UnitCell(...)

Define a new Bravais lattice type `MyBravaisLattice` with a unit cell constructor `UnitCell(expr)`.
If the notation is `N -> UnitCell(expr)`, the unit cell constructor will be dependent on the dimensionality `N`.
otherwise, the dimensionality will be inferred from the unit cell.
`N` is the dimensionality of the lattice.

## Examples

```jldoctest
julia> @bravaisdef MyBravaisLattice UnitCell([1 0; 0 1]);   # 2D square lattice

julia> MyBravaisLattice(3, 3)
9-site 2-dim MyBravaisLattice
boundaries: (none)
nnbonds:
1.0 =>
 1 => 1, [1, 0]
 1 => 1, [0, 1]
1.41421 =>
 1 => 1, [1, -1]
 1 => 1, [1, 1]
2.0 =>
 1 => 1, [2, 0]
 1 => 1, [0, 2]
```
"""
macro bravaisdef(type, expr)
    type_sym = type
    is_ndep = Meta.isexpr(expr, :->)
    if is_ndep
        @assert expr.args[1] isa Symbol "Invalid unitcell constructor; one parameter N expected"
        N_sym = expr.args[1]
        type_expr = :($type{$N_sym})
        type_def = :(const $type_expr = LatticeModels.BravaisLattice{$N_sym,
            <:LatticeModels.UnitCell{$(QuoteNode(type_sym))}})
        unitcell_construct = expr.args[2]
        unitcell_construct_signature =
            :(LatticeModels.construct_unitcell(::Type{$type_expr}) where $N_sym)
        lattice_constructor = quote
            function $type(sz::Vararg{LatticeModels.RangeT,$N_sym}; kw...) where $N_sym
                return LatticeModels.span_unitcells(
                    LatticeModels.construct_unitcell($type_expr), sz...; kw...)
            end
        end
    else
        type_expr = type
        type_def = :(const $type_expr = LatticeModels.BravaisLattice{N,
            <:LatticeModels.UnitCell{$(QuoteNode(type_sym))}} where N)
        unitcell_construct = expr
        unitcell_construct_signature = :(LatticeModels.construct_unitcell(::Type{$type_expr}))
        lattice_constructor = quote
            function $type(sz::Vararg{LatticeModels.RangeT}; kw...)
                return LatticeModels.span_unitcells(
                    LatticeModels.construct_unitcell($type_expr), sz...; kw...)
            end
        end
    end
    _add_parameter_to_call!(unitcell_construct, type_sym)
    return quote
        Core.@__doc__ $type_def
        $(Expr(:function, unitcell_construct_signature, unitcell_construct))
        $lattice_constructor
    end |> esc
end
function (::Type{T})(f::Function, args...; kw...) where T<:BravaisLattice
    filter(f, T(args...; kw...))
end

lattice_name(::Type{<:BravaisLattice{N, <:UnitCell{Sym}} where N}) where Sym = Sym
function (::Type{T})(::Val{LatticeRecipe}, args...; kw...) where {T<:BravaisLattice,LatticeRecipe}
    throw(ArgumentError("Recipe `$LatticeRecipe` not supported for $(lattice_name(T))"))
end
(::Type{T})(sym::Symbol, args...; kw...) where T<:BravaisLattice =
    T(Val(sym), args...; kw...)

"""
    SquareLattice{N}
Represents a square lattice in `N` dimensions.

---
    SquareLattice(sz::Int...)

Construct a square lattice of size `sz`.
"""
@bravaisdef SquareLattice N -> UnitCell(SMatrix{N,N}(I))
LatticeModels.site_coords(b::UnitCell{:squarelattice,N,1}, lp::BravaisPointer{N}) where {N} =
    vec(b.basissites) + lp.latcoords

"""
    TriangularLattice
Represents a triangular lattice.
Lattice vectors: `[1, 0]` and `[0.5, √3/2]`.

---
    TriangularLattice(a, b)

Construct a triangular lattice of a×b spanned unit cells.
"""
@bravaisdef TriangularLattice UnitCell([1 0.5; 0 √3/2])

"""
    TriangularLattice(:hex, hex_size[; center, kw...])

Construct a hexagon-shaped sample of a triangular lattice.

## Arguments:
- `hex_size`: the size of the hexagon.

## Keyword arguments:
- `center`: the center of the hexagon (in unit cell coordinates). Default is `(0, 0)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function TriangularLattice(::Val{:hex}, hex_size, center=(0, 0); kw...)
    j1ind = -hex_size + 1 + center[1]:hex_size - 1 + center[1]
    j2ind = -hex_size + 1 + center[2]:hex_size - 1 + center[2]
    TriangularLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords
        return -hex_size + 1 + sum(center) ≤ j1 + j2 ≤ hex_size - 1 + sum(center)
    end
end

"""
    TriangularLattice(:triangle, triangle_size[, center; kw...])

Construct a triangle-shaped sample of a triangular lattice.

## Arguments:
- `triangle_size`: the size of the triangle.

## Keyword arguments:
- `center`: the center of the triangle (in unit cell coordinates). Default is `(1, 1)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function TriangularLattice(::Val{:triangle}, triangle_size, center=(1, 1); kw...)
    j1ind = center[1]:triangle_size + center[1] - 1
    j2ind = center[2]:triangle_size + center[2] - 1
    TriangularLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords
        return j1 + j2 ≤ triangle_size - 1 + sum(center)
    end
end

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

"""
    HoneycombLattice(:hex, hex_size[, center; kw...])

Construct a hexagon-shaped sample of a honeycomb lattice.

## Arguments:
- `hex_size`: the size of the hexagon.
- `center`: the center of the hexagon (in unit cell coordinates). Default is `(0, 0)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function HoneycombLattice(::Val{:hex}, hex_size, center=(0, 0); kw...)
    j1ind = -hex_size + center[1]:hex_size - 1 + center[1]
    j2ind = -hex_size + 1 + center[2]:hex_size + center[2]
    HoneycombLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2
        if -hex_size ≤ j_prj ≤ hex_size
            return !(site.basindex == 2 && j_prj == hex_size ||
                site.basindex == 1 && j_prj == -hex_size)
        else
            return false
        end
    end
end

"""
    HoneycombLattice(:triangle, triangle_size[, center; kw...])

Construct a triangle-shaped sample of a honeycomb lattice.

## Arguments:
- `triangle_size`: the size of the triangle (number of hexagonal plaquettes).
- `center`: the center of the triangle (in unit cell coordinates). Default is `(1, 1)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function HoneycombLattice(::Val{:triangle}, triangle_size, center=(1, 1); preserve_tails=true, kw...)
    j1ind = center[1]:triangle_size + 1 + center[1]
    j2ind = center[2]:triangle_size + 1 + center[2]
    HoneycombLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2
        if j_prj ≤ triangle_size + 1
            if !preserve_tails && site.basindex == 1
                j1 == 0 && j2 == 0 && return false
                j1 == 0 && j2 == triangle_size + 1 && return false
                j1 == triangle_size + 1 && j2 == 0 && return false
            end
            return !(site.basindex == 2 && j_prj == triangle_size + 1)
        else
            return false
        end
    end
end

"""
    HoneycombLattice(:ribbon, r1, r2[, center; kw...])

Construct a graphene ribbon sample with zigzag edges.
To get armchair edges, simply rotate the lattice by 90 degrees.

## Arguments:
- `len`: the length of the ribbon.
- `wid`: the width of the ribbon.
- `center`: the unit cell coordinates of the bottom-left corner of the ribbon. Default is `(0, 0)`.

All other keyword arguments are passed to `span_unitcells` (see its documentation for details).
"""
function HoneycombLattice(::Val{:ribbon}, len, wid, center=(0, 0); kw...)
    j1ind = center[1] - wid:len + center[1] - 1
    j2ind = center[2]:wid + center[2] - 1
    HoneycombLattice(j1ind, j2ind; kw...) do site
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

"""
    KagomeLattice(:hex, hex_size[, center; kw...])

Construct a hexagon-shaped sample of a kagome lattice.

## Arguments:
- `hex_size`: the size of the hexagon.
- `center`: the center of the hexagon (in unit cell coordinates). Default is `(0, 0)`.
"""
function KagomeLattice(::Val{:hex}, hex_size, center=(0, 0); kw...)
    j1ind = -hex_size + 1 + center[1]:hex_size + center[1]
    j2ind = -hex_size + 1 + center[2]:hex_size + center[2]
    KagomeLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2 - 1
        if -hex_size ≤ j_prj ≤ hex_size - 1
            return !(site.basindex == 1 && j_prj == -hex_size ||
                    site.basindex == 2 && j1 == hex_size||
                    site.basindex == 3 && j2 == hex_size)
        else
            return false
        end
    end
end

"""
    KagomeLattice(:triangle, triangle_size[, center; kw...])

Construct a triangle-shaped sample of a kagome lattice.

## Arguments:
- `triangle_size`: the size of the triangle.
- `center`: the center of the triangle (in unit cell coordinates). Default is `(1, 1)`.
"""
function KagomeLattice(::Val{:triangle}, triangle_size, center=(1, 1); kw...)
    j1ind = center[1]:triangle_size + center[1]
    j2ind = center[2]:triangle_size + center[2]
    KagomeLattice(j1ind, j2ind; kw...) do site
        j1, j2 = site.latcoords .- center
        j_prj = j1 + j2
        if j_prj ≤ triangle_size
            return !(site.basindex == 1 && j_prj == 0 ||
                    site.basindex == 2 && j1 == triangle_size && j2 == 0 ||
                    site.basindex == 3 && j1 == 0 && j2 == triangle_size)
        else
            return false
        end
    end
end
