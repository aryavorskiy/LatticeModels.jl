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
    is_ndep = Meta.isexpr(expr, :->)
    if is_ndep
        @assert expr.args[1] isa Symbol "Invalid unitcell constructor; one parameter N expected"
        N_sym = expr.args[1]
        type_expr = :($type{$N_sym})
        type_def = :(const $type_expr = LatticeModels.BravaisLattice{$N_sym,
            <:LatticeModels.UnitCell{$(QuoteNode(type))}})
        unitcell_construct = expr.args[2]
        unitcell_construct_signature =
            :(LatticeModels.construct_unitcell(::Type{$type_expr}) where $N_sym)
        constructor_and_show = quote
            LatticeModels.construct_unitcell(::Type{$type}) =
                throw(ArgumentError("Unknown dimensionality for " * $(string(type)) *
                    "; Try using " * $(string(type)) * "{N}"))
            function $type(f::Function, sz::Vararg{LatticeModels.RangeT,$N_sym}; kw...) where $N_sym
                return LatticeModels.span_unitcells(
                    f, LatticeModels.construct_unitcell($type_expr), sz...; kw...)
            end
            Base.show(io::IO, ::Type{<:$type_expr}) where $N_sym =
                print(io, $(string(type)), "{", $N_sym, "}")
        end
    else
        type_def = :(const $type = LatticeModels.BravaisLattice{N,
            <:LatticeModels.UnitCell{$(QuoteNode(type))}} where N)
        unitcell_construct = expr
        unitcell_construct_signature = :(LatticeModels.construct_unitcell(::Type{$type}))
        constructor_and_show = quote
            function $type(f::Function, sz::Vararg{LatticeModels.RangeT}; kw...)
                return LatticeModels.span_unitcells(
                    f, LatticeModels.construct_unitcell($type), sz...; kw...)
            end
            Base.show(io::IO, ::Type{<:$type}) = print(io, $(string(type)))
        end
    end
    _add_parameter_to_call!(unitcell_construct, type)
    return quote
        Core.@__doc__ $type_def
        $(Expr(:function, unitcell_construct_signature, unitcell_construct))
        $constructor_and_show
        function Base.show(io::IO, ::MIME"text/plain", T::Type{<:$type})
            show(io, T)
            print(io, " (actually, ")
            Base._show_type(io, Base.inferencebarrier(T))
            print(io, ")")
        end
    end |> esc
end
function (::Type{T})(args...; kw...) where T<:BravaisLattice
    T(alwaystrue, T(args...; kw...))
end

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
    KagomeLattice
Represents a kagome lattice.

Lattice vectors: `[1, 0]` and `[0.5, √3/2]`,
three sites at `[0, 0]`, `[0.5, 0]` and `[0.25, √3/4]` in each unit cell.

---
    KagomeLattice(a, b)

Construct a kagome lattice of a×b spanned unit cells.
"""
@bravaisdef KagomeLattice UnitCell([1 0.5; 0 √3/2], [0 0.5 0.25; 0 0 √3/4])
