using SparseArrays
import QuantumOpticsBase

sites(l::AbstractLattice) = l
sites(l::BravaisLattice) = add_boundaries(l, BoundaryConditions())
struct LatticeBasis{LT<:AbstractLattice} <: QuantumOpticsBase.Basis
    shape::Int
    latt::LT
    function LatticeBasis(l::LT) where LT<:AbstractLattice
        new_l = sites(l)
        return new{typeof(new_l)}(length(l), new_l)
    end
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::AbstractSite) =
    basisstate(T, b, site_index(b.latt, site))
QuantumOpticsBase.basisstate(T::Type, l::AbstractLattice, site::AbstractSite) =
    basisstate(T, LatticeBasis(l), site)
QuantumOpticsBase.basisstate(l::AbstractLattice, site::AbstractSite) =
    basisstate(ComplexF64, l, site)

const CompositeLatticeBasis{S, BT, LT} = CompositeBasis{S, Tuple{BT, LatticeBasis{LT}}}
const AbstractLatticeBasis = Union{LatticeBasis, CompositeLatticeBasis}

const LatticeOperator = DataOperator{BT, BT} where BT<:LatticeBasis
const CompositeLatticeOperator = DataOperator{BT, BT} where BT<:CompositeLatticeBasis
const AbstractLatticeOperator = DataOperator{BT, BT} where BT<:AbstractLatticeBasis
