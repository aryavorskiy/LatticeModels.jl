"""
    SquareLattice{N}
Type alias for `Lattice{:square,N,1}`.

---
    SquareLattice(sz::Int...)

Constructs a square lattice of size `sz`.
"""
const SquareLattice{N} = BravaisLattice{N, <:UnitCell{:square,N,1}}
UnitCell{:square,N,1}() where N = UnitCell{:square}(SMatrix{N,N}(I))
default_bonds(::SquareLattice{N}, ::Val{1}) where {N} = Tuple(BravaisShift(axis=i) for i in 1:N)
default_bonds(::SquareLattice{N}, ::Val{2}) where {N} = Tuple(BravaisShift(one_hot(i, Val(N)) + k * one_hot(j, Val(N))) for i in 1:N for j in 1:i-1 for k in (-1, 1))
default_bonds(::SquareLattice{N}, ::Val{3}) where {N} = Tuple(BravaisShift(axis=i, dist=2) for i in 1:N)
LatticeModels.site_coords(b::UnitCell{:square,N,1}, lp::BravaisPointer{N}) where {N} =
    vec(b.basis) + lp.unit_cell

const TriangularLattice = BravaisLattice{2, <:UnitCell{:triangular,2,1}}
UnitCell{:triangular,2,1}() = UnitCell{:triangular}([1 0.5; 0 √3/2])
default_bonds(::TriangularLattice, ::Val{1}) = BravaisShift([0, 1]), BravaisShift([-1, 0]), BravaisShift([1, -1])
default_bonds(::TriangularLattice, ::Val{2}) = BravaisShift([1, 1]), BravaisShift([-2, 1]), BravaisShift([1, -2])
default_bonds(::TriangularLattice, ::Val{3}) = BravaisShift([0, 2]), BravaisShift([-2, 0]), BravaisShift([2, -2])

"""
    HoneycombLattice
Type alias for `Lattice{:honeycomb,2,2}`.

---
    HoneycombLattice(sz::Vararg{Int, 2})

Constructs a honeycomb lattice with a `sz`-size macrocell.
"""
const HoneycombLattice = BravaisLattice{2, <:UnitCell{:honeycomb,2,2}}
UnitCell{:honeycomb,2,2}() = UnitCell{:honeycomb}([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
default_bonds(::HoneycombLattice, ::Val{1}) = (BravaisShift(2 => 1), BravaisShift(2 => 1, axis=1), BravaisShift(2 => 1, axis=2))
default_bonds(::HoneycombLattice, ::Val{2}) = (
    BravaisShift(1 => 1, axis = 1),
    BravaisShift(2 => 2, axis = 1, dist=-1),
    BravaisShift(1 => 1, axis = 2, dist=-1),
    BravaisShift(2 => 2, axis = 2),
    BravaisShift(1 => 1, [-1, 1]),
    BravaisShift(2 => 2, [1, -1]))
default_bonds(::HoneycombLattice, ::Val{3}) = (
    BravaisShift(2 => 1, SA[1, 1]),
    BravaisShift(2 => 1, SA[1, -1]),
    BravaisShift(2 => 1, SA[-1, 1]))
