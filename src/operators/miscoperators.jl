import QuantumOpticsBase

function QuantumOpticsBase.diagonaloperator(bas::OneParticleBasis, lv::LatticeValue)
    check_samesites(lv, bas)
    N = internal_length(bas)
    return diagonaloperator(bas, repeat(lv.values, inner=N))
end
QuantumOpticsBase.diagonaloperator(lv::LatticeValue) =
    diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
QuantumOpticsBase.diagonaloperator(sys::OneParticleBasisSystem, lv::LatticeValue) =
    diagonaloperator(basis(sys), lv)
QuantumOpticsBase.diagonaloperator(bas::OneParticleBasis, crd::Union{SiteProperty,Symbol}) =
    diagonaloperator(bas, LatticeValue(lattice(bas), crd))
QuantumOpticsBase.diagonaloperator(sample::OneParticleBasisSystem, crd::Union{SiteProperty,Symbol}) =
    diagonaloperator(basis(sample), crd)
@accepts_system QuantumOpticsBase.diagonaloperator

"""
    coordoperators(sys)
    coordoperators(basis)
    coordoperators(lat[, internal])

Generate a `Tuple` of coordinate operators for the given lattice.

## Arguments
- `sys`: a `System` for which the coordinate operators are to be generated.
- `basis`: a one-particle `Basis` for which the coordinate operators are to be generated.
- `lat`: a lattice for which the coordinate operators are to be generated.
- `internal`: The basis for the internal degrees of freedom.
"""
coordoperators(lb::OneParticleBasis) =
    Tuple(diagonaloperator(lb, lv) for lv in coordvalues(lattice(lb)))
coordoperators(sample::OneParticleBasisSystem) = coordoperators(basis(sample))
@accepts_system coordoperators

"""
    coordoperator(sys, crd)
    coordoperator(basis, crd)
    coordoperator(lat[, internal], crd)

Generate a coordinate operator for the given lattice.

## Arguments
- `sys`: a `System` for which the coordinate operators are to be generated.
- `basis`: a one-particle `Basis` for which the coordinate operators are to be generated.
- `lat`: a lattice for which the coordinate operators are to be generated.
- `internal`: The basis for the internal degrees of freedom.
- `crd`: The coordinate to generate the operator for. Must be an integer representing the
    coordinate index (e. g. `1` for x, `2` for y, etc.) or a symbol (e. g. `:x`, `:y`, etc.).
"""
coordoperator(lb::OneParticleBasis, crd) = coordoperator(sample(lb), crd)
coordoperator(sample::Sample, i::Int) = diagonaloperator(sample, Coord(i))
coordoperator(sample::Sample, sym::Symbol) = diagonaloperator(sample, SitePropertyAlias{sym}())
@accepts_system coordoperator

const SampleIndex = Union{AbstractSite, Tuple{AbstractSite, Int}}
function to_index(sample::Sample, site::AbstractSite)
    hasinternal(sample) &&
        throw(ArgumentError("Basis has internal degrees of freedom, specify the internal index"))
    return site_index(lattice(sample), site)
end
function to_index(sample::Sample, ind::Tuple{AbstractSite, Int})
    site, i = ind
    N = internal_length(sample)
    i > N && throw(ArgumentError("Index $i is out of bounds for the internal basis of length $N"))
    return (site_index(lattice(sample), site) - 1) * N + i
end
to_index(any, ind::SampleIndex) = to_index(sample(any), ind)

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, ind::SampleIndex) =
    basisstate(T, b, to_index(b, ind))
QuantumOpticsBase.basisstate(T::Type, op::OneParticleBasisSystem, ind::SampleIndex) =
    basisstate(T, basis(op), ind)

function QuantumOpticsBase.transition(sys::System, site1::AbstractSite, site2::AbstractSite)
    builder = FastOperatorBuilder(sys)
    builder[site1, site2] = internal_one(sys)
    return Operator(builder)
end
function QuantumOpticsBase.transition(b::Basis, ind1::SampleIndex, ind2::SampleIndex)
    return transition(b, to_index(b, ind1), to_index(b, ind2))
end
function QuantumOpticsBase.transition(sys::System, ind1::SampleIndex, ind2::SampleIndex)
    return transition(basis(sys), ind1, ind2)
end
@accepts_system QuantumOpticsBase.transition

QuantumOpticsBase.number(b::Basis, ind1::SampleIndex) = QuantumOpticsBase.transition(b, ind1, ind1)
QuantumOpticsBase.number(sys::System, ind1::SampleIndex) = QuantumOpticsBase.transition(sys, ind1, ind1)
@accepts_system QuantumOpticsBase.number

QuantumOpticsBase.create(b::Basis, ind1::SampleIndex) = create(b, to_index(b, ind1))
QuantumOpticsBase.create(sys::System, ind1::SampleIndex) = create(basis(sys), ind1)
@accepts_system QuantumOpticsBase.create

QuantumOpticsBase.destroy(b::Basis, ind1::SampleIndex) = destroy(b, to_index(b, ind1))
QuantumOpticsBase.destroy(sys::System, ind1::SampleIndex) = destroy(basis(sys), ind1)
@accepts_system QuantumOpticsBase.destroy
