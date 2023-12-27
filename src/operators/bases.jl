using SparseArrays
import QuantumOpticsBase
import QuantumOpticsBase: basis

struct LatticeBasis{ST<:Sites} <: QuantumOpticsBase.Basis
    shape::Int
    sites::ST
    function LatticeBasis(l::LT) where LT<:AbstractLattice
        new_l = sites(l)
        return new{typeof(new_l)}(length(l), new_l)
    end
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.sites == lb2.sites

function Base.show(io::IO, ::MIME"text/plain", bas::LatticeBasis)
    print(io, "LatticeBasis(")
    summary(io, bas.sites.latt)
    print(io, ")")
end

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::AbstractSite) =
    basisstate(T, b, site_index(b.sites, site))
QuantumOpticsBase.basisstate(T::Type, l::AbstractLattice, site::AbstractSite) =
    basisstate(T, LatticeBasis(l), site)
QuantumOpticsBase.basisstate(l::AbstractLattice, site::AbstractSite) =
    basisstate(ComplexF64, l, site)

const CompositeLatticeBasis{S, BT, LT} = CompositeBasis{S, Tuple{BT, LatticeBasis{LT}}}
const AbstractLatticeBasis = Union{LatticeBasis, CompositeLatticeBasis}

const LatticeOperator = DataOperator{BT, BT} where BT<:LatticeBasis
const CompositeLatticeOperator = DataOperator{BT, BT} where BT<:CompositeLatticeBasis
const AbstractLatticeOperator = DataOperator{BT, BT} where BT<:AbstractLatticeBasis

sites(lb::LatticeBasis) = lb.sites
sites(lb::CompositeLatticeBasis) = sites(lb.bases[end])
sites(op::AbstractOperator) = sites(basis(op))
sample(op::AbstractOperator) = sample(basis(op))
