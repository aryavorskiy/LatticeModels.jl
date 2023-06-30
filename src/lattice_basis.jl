import Base: ==
import QuantumOpticsBase: diagonaloperator

struct LatticeBasis{LT<:Lattice} <: QuantumOpticsBase.Basis
    shape::SVector{Int, 1}
    latt::LT
    LatticeBasis(l::LT) where LT<:Lattice = new{LT}(SA[length(l)], l)
end
==(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt
lattice(lb::LatticeBasis) = lb.latt

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::LatticeSite) =
    QuantumOpticsBase.basisstate(T, b, site_index(b.latt, site))
function QuantumOpticsBase.diagonaloperator(lv::LatticeValue)
    QuantumOpticsBase.diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
end
function QuantumOpticsBase.diagonaloperator(f::Function, b::LatticeBasis)
    QuantumOpticsBase.diagonaloperator(b, f.(b.latt))
end

"""
coord_operators(lb::LatticeBasis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
coord_operators(lb::LatticeBasis) = (diagonaloperator(lv) for lv in coord_values(lb.latt))
