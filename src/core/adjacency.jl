"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`AdjacencyMatrix`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `AdjacencyMatrix` which connects sites that were connected by `≤n` bonds of the previous `AdjacencyMatrix`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct AdjacencyMatrix{LT<:BravaisLattice}
    lattice::LT
    bmat::Matrix{Bool}
    function AdjacencyMatrix(l::LT, bmat) where {LT<:BravaisLattice}
        !all(size(bmat) .== length(l)) && error("inconsistent connectivity matrix size")
        new{LT}(l, bmat .| bmat' .| Matrix(I, length(l), length(l)))
    end
    function AdjacencyMatrix(l::BravaisLattice)
        AdjacencyMatrix(l, Matrix(I, length(l), length(l)))
    end
end

lattice(bs::AdjacencyMatrix) = bs.lattice
match(bs::AdjacencyMatrix, site1::BravaisSite, site2::BravaisSite) =
    bs.bmat[site_index(lattice(bs), site1), site_index(lattice(bs), site2)]

function Base.:(|)(bss1::AdjacencyMatrix, bss2::AdjacencyMatrix)
    check_samelattice(bss1, bss2)
    AdjacencyMatrix(lattice(bss1), bss1.bmat .| bss2.bmat)
end

Base.:(^)(bs1::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(bs1.lattice, bs1.bmat^n .!= 0)
