QuantumOpticsBase.diagonaloperator(lv::LatticeValue) =
    QuantumOpticsBase.diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
function QuantumOpticsBase.diagonaloperator(lb::AbstractLatticeBasis, lv::LatticeValue)
    check_samelattice(lv, lb)
    N = internal_length(lb)
    return diagonaloperator(lb, repeat(lv.values, inner=N))
end
QuantumOpticsBase.diagonaloperator(sys::OneParticleBasisSystem, lv::LatticeValue) =
    QuantumOpticsBase.diagonaloperator(basis(sys), lv)
@accepts_system QuantumOpticsBase.diagonaloperator

"""
    coord_operators(l::Lattice[, ib::Basis])
    coord_operators(lb::AbstractLatticeBasis)

Generate a `Tuple` of coordinate operators for given `sample`.

Standard rules for functions accepting `Sample`s apply.
"""
coord_operators(lb::AbstractLatticeBasis) =
    Tuple(diagonaloperator(lb, lv) for lv in coord_values(lattice(lb)))
coord_operators(sample::OneParticleBasisSystem) = coord_operators(basis(sample))
@accepts_system coord_operators

param_operator(lb::AbstractLatticeBasis, crd) =
    diagonaloperator(lb, param_value(lattice(lb), crd))
param_operator(sample::OneParticleBasisSystem, crd) = param_operator(basis(sample), crd)
@accepts_system param_operator

"""
    QuantumOpticsBase.transition(sys::System, site1::LatticeSite, site2::LatticeSite[, op; field])
    QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int[, op; field])

Generate a transition operator between two local states in lattice space.
States can be defined by `LatticeSite`s or integers.

Standard rules for functions accepting `System`s apply.
"""
function QuantumOpticsBase.transition(sys::System, site1::AbstractSite, site2::AbstractSite, op=internal_one(sample); field=NoField())
    return build_operator(sys, op => site1 => site2, field=field)
end
QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int, op=internal_one(sys); field=NoField()) =
    build_operator(sys, op => lattice(sys)[i1] => lattice(sample)[i2], field=field)
@accepts_system QuantumOpticsBase.transition
