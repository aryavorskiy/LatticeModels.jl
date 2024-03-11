using SparseArrays
import QuantumOpticsBase
import QuantumOpticsBase: basis

"""
    LatticeBasis <: QuantumOpticsBase.Basis

Basis for a lattice system.
"""
struct LatticeBasis{LT} <: QuantumOpticsBase.Basis
    shape::Int
    lat::LT
    function LatticeBasis(l::LT) where LT<:AbstractLattice
        new_l = stripparams(l)
        hasparam(l, :latticetype) && (new_l = settype(new_l, gettype(l)))
        return new{typeof(new_l)}(length(l), new_l)
    end
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.lat == lb2.lat

function Base.show(io::IO, ::MIME"text/plain", bas::LatticeBasis)
    print(io, "LatticeBasis(")
    summary(io, bas.lat)
    print(io, ")")
end

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::AbstractSite) =
    basisstate(T, b, site_index(b.lat, site))
QuantumOpticsBase.basisstate(T::Type, l::AbstractLattice, site::AbstractSite) =
    basisstate(T, LatticeBasis(l), site)
QuantumOpticsBase.basisstate(l::AbstractLattice, site::AbstractSite) =
    basisstate(ComplexF64, l, site)

const CompositeLatticeBasis{S, BT, LT} = CompositeBasis{S, Tuple{BT, LatticeBasis{LT}}}
const OneParticleBasis = Union{LatticeBasis, CompositeLatticeBasis}
const AbstractLatticeBasis = Union{OneParticleBasis, ManyBodyBasis{<:OneParticleBasis}}

const LatticeOperator = DataOperator{BT, BT} where BT<:LatticeBasis
const CompositeLatticeOperator = DataOperator{BT, BT} where BT<:CompositeLatticeBasis
const OneParticleOperator = DataOperator{BT, BT} where BT<:OneParticleBasis
const AbstractLatticeOperator = DataOperator{BT, BT} where BT<:AbstractLatticeBasis

const StateType{BT<:Basis} = Union{DataOperator{BT,BT}, Ket{BT}, Bra{BT}}
# These function interpret any state type (vector/matrix) as a density matrix
@inline matrix_element(op::DataOperator, i, j) = op.data[i, j]
@inline matrix_element(ket::Ket, i, j) = ket.data[i] * ket.data[j]'
@inline matrix_element(bra::Bra, i, j) = bra.data[j] * bra.data[i]'

lattice(lb::LatticeBasis) = lb.lat
lattice(lb::CompositeLatticeBasis) = lattice(lb.bases[end])
lattice(ms::ManyBodyBasis{<:OneParticleBasis}) = lattice(ms.onebodybasis)
lattice(state::StateType) = lattice(basis(state))

LatticeValue(ket::Ket{<:LatticeBasis}) = LatticeValue(lattice(ket), ket.data)
LatticeValue(bra::Bra{<:LatticeBasis}) = LatticeValue(lattice(bra), bra.data)
LatticeValue(::StateVector{<:CompositeLatticeBasis}) =
    throw(ArgumentError("""Cannot convert a state on a lattice basis \
with on-site degrees of freedom to a lattice value"""))
