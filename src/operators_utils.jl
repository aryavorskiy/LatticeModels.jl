import QuantumOpticsBase: basis, samebases, check_samebases

function QuantumOpticsBase.diagonaloperator(lv::LatticeValue)
    diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
end
function QuantumOpticsBase.diagonaloperator(f::Function, b::LatticeBasis)
    diagonaloperator(b, f.(b.latt))
end

"""
coord_operators(lb::LatticeBasis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
coord_operators(lb::LatticeBasis) = Tuple(diagonaloperator(lb, lv.values) for lv in coord_values(lb.latt))
function coord_operators(cb::CompositeLatticeBasis)
    N = length(cb.bases[1])
    lb = cb.bases[2]
    return Tuple(diagonaloperator(cb, repeat(lv.values, inner=N)) for lv in coord_values(lb.latt))
end
coord_operators(l::Lattice) = coord_operators(LatticeBasis(l))
coord(lb::LatticeBasis, crd) = diagonaloperator(lb, [getproperty(site, crd) for site in lb.latt])
coord(l::Lattice, crd) = coord(LatticeBasis(l), crd)

function site_density(ket::Ket{<:LatticeBasis})
    LatticeValue(lattice(ket), map(abs2, ket.data))
end

function site_density(ket::Ket{<:CompositeLatticeBasis})
    l = lattice(ket)
    N = length(internal_basis(ket))
    LatticeValue(l, [sum(abs2, @view(ket.data[(i - 1) * N + 1: i * N])) for i in 1:length(l)])
end

site_density(bra::Bra) = site_density(dagger(bra))

function site_density(op::LatticeOperator)
    LatticeValue(lattice(op), diag(op.data))
end

function site_density(op::CompositeLatticeOperator)
    l = lattice(op)
    N = length(internal_basis(op))
    dg = diag(op.data)
    LatticeValue(l, [sum(@view(dg[(i - 1) * N + 1: i * N])) for i in 1:length(l)])
end
