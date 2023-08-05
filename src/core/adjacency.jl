"""
    AbstractGraph <: Function

A function-like object that accepts a `Lattice` and two `LatticeSite`s and returns whether the site pair is *selected* or not.
Implements a built-in sanity check algorithm to make sure the pair set was defined on a correct lattice.

Define these functions for all subclasses:
- `LatticeModels.lattice(::YourSelector)` must return the lattice your selector was defined on.
- `LatticeModels.match(::YourSelector, site1::LatticeSite, site2::LatticeSite)` must return whether the `(site1, site2)` pair is *selected*.
"""
abstract type AbstractGraph end
match(ps::AbstractGraph, ::LatticeSite, ::LatticeSite) =
    error("match(::$(typeof(ps)), ::LatticeSite, ::LatticeSite) must be explicitly implemented")
lattice(ps::AbstractGraph) = error("lattice(::$(typeof(ps))) must be explicitly implemented")

check_lattice_fits(::Any, ::Lattice) = nothing
check_lattice_fits(ps::AbstractGraph, l::Lattice) = check_issublattice(l, lattice(ps))

"""
    Domains(domains::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `domains` is the same on two sites.
"""
struct Domains{LT} <: AbstractGraph
    domains::LT
    Domains(domains::LT) where LT<:LatticeValue = new{LT}(domains)
end
lattice(ps::Domains) = lattice(ps.domains)
match(ps::Domains, site1::LatticeSite, site2::LatticeSite) =
    ps.domains[site1] == ps.domains[site2]

"""
    PairLhsGraph(lhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `lhs` is true on the first site of the pair.
"""
struct PairLhsGraph <: AbstractGraph
    lhs::LatticeValue
end
lattice(ps::PairLhsGraph) = lattice(ps.lhs)
match(ps::PairLhsGraph, site1::LatticeSite, ::LatticeSite) = ps.lhs[site1]

"""
    PairRhsGraph(rhs::LatticeValue)

A selector used for hopping operator definition or currents materialization.

Takes a `LatticeValue`.
A pair matches the selector if the value of `rhs` is true on the first site of the pair.
"""
struct PairRhsGraph <: AbstractGraph
    rhs::LatticeValue
end
lattice(ps::PairRhsGraph) = lattice(ps.rhs)
match(ps::PairRhsGraph, ::LatticeSite, site2::LatticeSite) = ps.rhs[site2]

"""
    AdjacencyMatrix{LT} where {LT<:Lattice}

Represents the bonds on some lattice.

`AdjacencyMatrix`s can be combined with the `|` operator and negated with the `!` operator.
Also you can create a `AdjacencyMatrix` which connects sites that were connected by `≤n` bonds of the previous `AdjacencyMatrix`
by taking its power: `bs2 = bs1 ^ n`.
"""
struct AdjacencyMatrix{LT<:Lattice} <: AbstractGraph
    lattice::LT
    bmat::Matrix{Bool}
    function AdjacencyMatrix(l::LT, bmat) where {LT<:Lattice}
        !all(size(bmat) .== length(l)) && error("inconsistent connectivity matrix size")
        new{LT}(l, bmat .| bmat' .| Matrix(I, length(l), length(l)))
    end
    function AdjacencyMatrix(l::Lattice)
        AdjacencyMatrix(l, Matrix(I, length(l), length(l)))
    end
end

lattice(bs::AdjacencyMatrix) = bs.lattice
match(bs::AdjacencyMatrix, site1::LatticeSite, site2::LatticeSite) =
    bs.bmat[site_index(lattice(bs), site1), site_index(lattice(bs), site2)]

function Base.:(|)(bss1::AdjacencyMatrix, bss2::AdjacencyMatrix)
    check_samelattice(bss1, bss2)
    AdjacencyMatrix(lattice(bss1), bss1.bmat .| bss2.bmat)
end

Base.:(^)(bs1::AdjacencyMatrix, n::Int) =
    AdjacencyMatrix(bs1.lattice, bs1.bmat^n .!= 0)

struct InvertedGraph{GT} <: AbstractGraph
    graph::GT
end
lattice(ig::InvertedGraph) = lattice(ig.graph)
match(ig::InvertedGraph, site1, site2) = match(ig.graph, site1, site2)
Base.:(!)(g::AbstractGraph) = InvertedGraph(g)
Base.:(!)(g::InvertedGraph) = g.graph
Base.:(!)(bs::AdjacencyMatrix) =
    AdjacencyMatrix(bs.lattice, Matrix(I, length(bs.lattice), length(bs.lattice)) .| .!(bs.bmat))
