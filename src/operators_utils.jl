import QuantumOpticsBase: basis, samebases, check_samebases

QuantumOpticsBase.diagonaloperator(lv::LatticeValue) =
    QuantumOpticsBase.diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
function QuantumOpticsBase.diagonaloperator(lb::AbstractLatticeBasis, lv::LatticeValue)
    check_samelattice(lv, lb)
    N = internal_length(lb)
    return diagonaloperator(lb, repeat(lv.values, inner=N))
end
QuantumOpticsBase.diagonaloperator(sample::Sample, lv::LatticeValue) =
    QuantumOpticsBase.diagonaloperator(basis(sample), lv)
@accepts_sample QuantumOpticsBase.diagonaloperator

"""
    coord_operators(sample::Sample)
    coord_operators(l::Lattice[, ib::Basis])
    coord_operators(lb::AbstractLatticeBasis)

Generate a `Tuple` of coordinate operators for given `sample`.

Standard rules for functions accepting `Sample`s apply.
"""
coord_operators(lb::AbstractLatticeBasis) =
    Tuple(diagonaloperator(lb, lv) for lv in coord_values(lattice(lb)))
coord_operators(sample::Sample) = coord_operators(basis(sample))
@accepts_sample coord_operators

coord_operator(lb::AbstractLatticeBasis, crd) =
    diagonaloperator(lb, coord_value(lattice(lb), crd))
coord_operator(sample::Sample, crd) = coord_operator(basis(sample), crd)
@accepts_sample coord_operator

"""
    QuantumOpticsBase.transition(sys::System, site1::LatticeSite, site2::LatticeSite[, op; field])
    QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int[, op; field])

Generate a transition operator between two local states in lattice space.
States can be defined by `LatticeSite`s or integers.

Standard rules for functions accepting `System`s apply.
"""
function QuantumOpticsBase.transition(sys::System, site1::BravaisSite, site2::BravaisSite, op=internal_one(sample); field=NoField())
    return build_operator(sys, op => site1 => site2, field=field)
end
QuantumOpticsBase.transition(sys::System, i1::Int, i2::Int, op=internal_one(sample); field=NoField()) =
    build_operator(sys, op => lattice(sample)[i1] => lattice(sample)[i2], field=field)
@accepts_system QuantumOpticsBase.transition

"""
    lattice_density(ket::Ket)
    lattice_density(bra::Bra)
    lattice_density(densitymat::Operator)

Calculate local density of given state.
The state can be expressed as a ket/bra vector or a density matrix.
The output of this function is a `LatticeValue`.
"""
function lattice_density(ket::Ket{<:LatticeBasis})
    LatticeValue(lattice(ket), map(abs2, ket.data))
end
function lattice_density(ket::Ket{<:CompositeLatticeBasis})
    l = lattice(ket)
    N = internal_length(ket)
    LatticeValue(l, [sum(abs2, @view(ket.data[(i - 1) * N + 1: i * N])) for i in 1:length(l)])
end

lattice_density(bra::Bra) = lattice_density(dagger(bra))

function lattice_density(op::LatticeOperator)
    LatticeValue(lattice(op), real.(diag(op.data)))
end

function lattice_density(op::CompositeLatticeOperator)
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

function QuantumOpticsBase.ptrace(op::LatticeModels.CompositeLatticeOperator, sym::Symbol)
    bs = length(basis(op).bases)
    if sym === :lattice
        return ptrace(op, bs)
    elseif sym === :internal
        return ptrace(op, Tuple(1:bs-1))
    else
        throw(ArgumentError("Invalid subspace symbol ':$sym'; ':lattice' or ':internal' expected"))
    end
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

get_site_periodic(::BravaisLattice, ::Nothing) = nothing
function get_site_periodic(l::BravaisLattice, site::BravaisPointer)
    new_site = get_site(l, site)
    new_site === nothing ? nothing : shift_site(PeriodicBoundaryConditions(), l, new_site)[2]
end
function adjacency_matrix(l::BravaisLattice, bss::SiteOffset...)
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
