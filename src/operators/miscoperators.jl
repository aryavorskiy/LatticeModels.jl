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

coordoperator(lb::OneParticleBasis, crd) = coordoperator(sample(lb), crd)
coordoperator(sample::Sample, i::Int) = diagonaloperator(sample, Coord(i))
coordoperator(sample::Sample, sym::Symbol) = diagonaloperator(sample, SitePropertyAlias{sym}())
@accepts_system coordoperator

"""
    QuantumOpticsBase.transition(sys::System, site1::LatticeSite, site2::LatticeSite[, op; field])
    QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int[, op; field])

Generate a transition operator between two local states in lattice space.
States can be defined by `LatticeSite`s or integers.

Standard rules for functions accepting `System`s apply.
"""
function QuantumOpticsBase.transition(sys::System, site1::AbstractSite, site2::AbstractSite, op=internal_one(sample); field=NoField())
    return construct_operator(sys, op => site1 => site2, field=field)
end
QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int, op=internal_one(sys); field=NoField()) =
    construct_operator(sys, op => lattice(sys)[i1] => lattice(sample)[i2], field=field)
@accepts_system QuantumOpticsBase.transition
