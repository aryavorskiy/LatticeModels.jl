using SparseArrays
import QuantumOpticsBase

struct LatticeBasis{LT<:BravaisLattice} <: QuantumOpticsBase.Basis
    shape::Int
    latt::LT
    LatticeBasis(l::LT) where LT<:BravaisLattice = new{LT}(length(l), l)
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::BravaisSite) =
    basisstate(T, b, site_index(b.latt, site))

const CompositeLatticeBasis{S, BT, LT} = CompositeBasis{S, Tuple{BT, LatticeBasis{LT}}}
const AbstractLatticeBasis = Union{LatticeBasis, CompositeLatticeBasis}

const LatticeOperator = DataOperator{BT, BT} where BT<:LatticeBasis
const CompositeLatticeOperator = DataOperator{BT, BT} where BT<:CompositeLatticeBasis
const AbstractLatticeOperator = DataOperator{BT, BT} where BT<:AbstractLatticeBasis
