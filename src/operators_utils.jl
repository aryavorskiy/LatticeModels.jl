import QuantumOpticsBase: basis, samebases, check_samebases

QuantumOpticsBase.diagonaloperator(lv::LatticeValue) =
    diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
function QuantumOpticsBase.diagonaloperator(lb::AbstractLatticeBasis, lv::LatticeValue)
    check_lattice_match(lv, lb)
    N = internal_length(lb)
    return diagonaloperator(lb, repeat(lv.values, inner=N))
end
QuantumOpticsBase.diagonaloperator(sample::Sample, lv::LatticeValue) =
    QuantumOpticsBase.diagonaloperator(basis(sample), lv)
@accepts_sample QuantumOpticsBase.diagonaloperator

"""
coord_operators(sample::Sample)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
coord_operators(lb::AbstractLatticeBasis) =
    Tuple(diagonaloperator(lb, lv) for lv in coord_values(lattice(lb)))
coord_operators(sample::Sample) = coord_operators(basis(sample))
@accepts_sample coord_operators

coord_operator(lb::AbstractLatticeBasis, crd) =
    diagonaloperator(lb, coord_value(lattice(lb), crd))
coord_operator(sample::Sample, crd) = coord_operator(basis(sample), crd)
@accepts_sample coord_operator

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
    LatticeValue(lattice(op), real.(diag(op.data)))
end

function site_density(op::CompositeLatticeOperator)
    l = lattice(op)
    N = internal_length(op)
    dg = diag(op.data)
    LatticeValue(l, [real(sum(@view(dg[(i - 1) * N + 1: i * N]))) for i in 1:length(l)])
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

function get_site_periodic(l::Lattice, site::LatticePointer)
    new_site = get_site(l, site)
    new_site === nothing ? nothing : shift_site(PeriodicBoundaryConditions(), l, new_site)[2]
end
function adjacency_matrix(l::Lattice, bss::SiteOffset...)
    matrix = zeros(Bool, length(l), length(l))
    for bs in bss
        for (i, site) in enumerate(l)
            new_site = get_site_periodic(l, site + bs)
            j = site_index(l, new_site)
            j === nothing && continue
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
