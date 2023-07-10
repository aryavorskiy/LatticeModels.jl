using StaticArrays

const SingleBond{N} = Pair{LatticeSite{N}, LatticeSite{N}}

"""
    Bonds{N}

A struct representing bonds in some direction in a lattice.
- `site_indices`: a `NTuple{2, Int}` with indices of sites connected by the bond.
- `translate_uc`: the unit cell offset.
"""
struct Bonds{N}
    site_indices::Pair{Int,Int}
    translate_uc::SVector{N, Int}
end

"""
    Bonds(site_indices; kwargs...)

A convenient constructor for a `Bonds` object.
`site_indices` is a `::Int => ::Int` pair with indices of sites connected by the bond; `1 => 1` is the default value.

**Keyword arguments:**
- `translate_uc`: the unit cell offset, a `Tuple` or a `SVector`. Zeros by default.
- `axis`: overrides `translate_uc` and sets its components to zero on all axes except given.

If `site_indices` are equal and `translate_uc` is zero, this means that the bond connects each site with itself,
in which case an error will be thrown.
Note that the dimension count for the bond is static, it is automatically compatible to higher-dimensional lattices.
"""
function Bonds(site_indices=Pair(1, 1); kw...)
    tr_vc = if :axis in keys(kw)
        one_hot(kw[:axis], kw[:axis])
    elseif :translate_uc in keys(kw)
        kw[:translate_uc]
    elseif site_indices[1] != site_indices[2]
        SA[0]
    end
    if iszero(tr_vc) && ==(site_indices...)
        throw(ArgumentError("bond connects site to itself"))
    end
    Bonds(site_indices, tr_vc)
end

Base.:(==)(h1::Bonds, h2::Bonds) =
    all(getproperty(h1, fn) == getproperty(h2, fn) for fn in fieldnames(Bonds))

function Base.show(io::IO, ::MIME"text/plain", hop::Bonds)
    println(io, "Bonds connect site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.translate_uc)")
end
dims(::Bonds{N}) where N = N

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
check_lattice_fits(ps::AbstractGraph, l::Lattice) = check_is_sublattice(l, lattice(ps))

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

struct InvertedGraph{GT} <: AbstractGraph
    graph::GT
end
lattice(ig::InvertedGraph) = lattice(ig.graph)
match(ig::InvertedGraph, site1, site2) = match(ig.graph, site1, site2)
Base.:(!)(g::AbstractGraph) = InvertedGraph(g)
Base.:(!)(g::InvertedGraph) = g.graph
