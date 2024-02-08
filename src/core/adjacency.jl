"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`AdjacencyMatrix`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `AdjacencyMatrix` which connects sites that were connected by `≤n` bonds of the previous `AdjacencyMatrix`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct AdjacencyMatrix{LT<:Sites}
    sites::LT
    bmat::Matrix{Bool}
    function AdjacencyMatrix(l::LT, connectivity_matrix) where {LT<:AbstractLattice}
        @check_size connectivity_matrix :square
        @check_size l size(connectivity_matrix, 1)
        s = sites(l)
        new{typeof(s)}(s,
            connectivity_matrix .| connectivity_matrix' .| Matrix(I, length(l), length(l)))
    end
    function AdjacencyMatrix(l::AbstractLattice)
        AdjacencyMatrix(l, Matrix(I, length(l), length(l)))
    end
end

sites(bs::AdjacencyMatrix) = bs.sites
Base.getindex(bs::AdjacencyMatrix, site1::AbstractSite, site2::AbstractSite) =
    bs.bmat[site_index(bs.sites, site1), site_index(bs.sites, site2)]

function Base.:(|)(bss1::AdjacencyMatrix, bss2::AdjacencyMatrix)
    check_samesites(bss1, bss2)
    AdjacencyMatrix(bss1.sites, bss1.bmat .| bss2.bmat)
end

Base.:(^)(bs1::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(bs1.sites, bs1.bmat^n .!= 0)

Base.:(!)(bs::AdjacencyMatrix) = AdjacencyMatrix(bs.sites, .!bs.bmat)
