import Base: ==
import QuantumOpticsBase

struct LatticeBasis{LT<:Lattice} <: QuantumOpticsBase.Basis
    shape::SVector{1, Int}
    latt::LT
    LatticeBasis(l::LT) where LT<:Lattice = new{LT}(SA[length(l)], l)
end
==(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt
lattice(lb::LatticeBasis) = lb.latt

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::LatticeSite) =
    basisstate(T, b, site_index(b.latt, site))
function QuantumOpticsBase.diagonaloperator(lv::LatticeValue)
    diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
end
function QuantumOpticsBase.diagonaloperator(f::Function, b::LatticeBasis)
    diagonaloperator(b, f.(b.latt))
end

function densityoperator(lb::LatticeBasis, l::Lattice)
    check_is_sublattice(lb.latt, l)
    diagonaloperator(in(l), lb)
end

"""
coord_operators(lb::LatticeBasis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
coord_operators(lb::LatticeBasis) = Tuple(diagonaloperator(lb, lv) for lv in coord_values(lb.latt))
# coord(lb::LatticeBasis, coord) = diagonaloperator(lb, coord_values(lb.latt)[_parse_axis_descriptor(coord)])

const LatticeOperator{T, MT} = Operator{T, T, MT} where T<:LatticeBasis
