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
    N = internal_length(ket)
    LatticeValue(l, [sum(abs2, @view(ket.data[(i - 1) * N + 1: i * N])) for i in 1:length(l)])
end

site_density(bra::Bra) = site_density(dagger(bra))

function site_density(op::LatticeOperator)
    LatticeValue(lattice(op), diag(op.data))
end

function site_density(op::CompositeLatticeOperator)
    l = lattice(op)
    N = internal_length(op)
    dg = diag(op.data)
    LatticeValue(l, [sum(@view(dg[(i - 1) * N + 1: i * N])) for i in 1:length(l)])
end

function diag_reduce(f, op::AbstractLatticeOperator)
    l = lattice(op)
    N = internal_length(op)
    LatticeValue(l,
        [f(@view op.data[(i - 1) * N + 1: i * N, (i - 1) * N + 1: i * N]
        ) for i in 1:length(l)])
end

"""
    adjacency_matrix(op::Operator)

Generates an `AdjacencyMatrix` for the provided operator.
"""
function adjacency_matrix(op::LatticeOperator)
    matrix = Bool[!iszero(op.data[i, j])
                  for i in 1:length(lattice(op)), j in 1:length(lattice(op))]
    return AdjacencyMatrix(lattice(op), matrix)
end
function adjacency_matrix(op::CompositeLatticeOperator)
    n = internal_length(op)
    ind(k) = (k - 1) * n + 1 : k * n
    matrix = Bool[!iszero(op.data[ind(i), ind(j)])
                  for i in 1:length(lattice(op)), j in 1:length(lattice(op))]
    return AdjacencyMatrix(lattice(op), matrix)
end

function adjacency_matrix(l::Lattice, bss::Bonds...)
    matrix = zeros(Bool, length(l), length(l))
    for bs in bss
        for (i, site) in enumerate(l)
            j = site_index(l, site + LatticeOffset(l, bs))
            matrix[i, j] = matrix[j, i] = true
        end
    end
    return AdjacencyMatrix(l, matrix)
end

function apply_field!(op::AbstractLatticeOperator, field::AbstractField)
    l = lattice(op)
    N = internal_length(op)
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            op.data[(i - 1) * N + 1: i * N, (j - 1) * N + 1: j * N] *=
                exp(-2π * im * line_integral(field, site1.coords, site2.coords))
        end
    end
end
