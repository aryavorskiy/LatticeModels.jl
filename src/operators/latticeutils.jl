import QuantumOpticsBase: basis, samebases, check_samebases

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

function lattice_density(ket::Ket{<:ManyBodyBasis{<:Any,<:AbstractLatticeBasis}})
    vs = zeros(length(basis(ket).onebodybasis))
    for i in 1:length(ket)
        occ = basis(ket).occupations[i]
        @. vs += occ * abs2(ket.data[i])
    end
    N = internal_length(ket)
    l = lattice(ket)
    LatticeValue(l, [@view(vs[(i - 1) * N + 1: i * N]) for i in 1:length(l)])
end

function lattice_density(op::DataOperator{BT, BT} where BT<:ManyBodyBasis{<:Any, <:AbstractLatticeBasis})
    vs = zeros(length(basis(op).onebodybasis))
    ds = diag(op.data)
    for i in eachindex(ds)
        occ = basis(op).occupations[i]
        @. vs += occ * ds[i]
    end
    N = internal_length(op)
    l = lattice(op)
    LatticeValue(l, [@view(vs[(i - 1) * N + 1: i * N]) for i in 1:length(l)])
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
