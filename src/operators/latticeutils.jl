import QuantumOpticsBase: basis, samebases, check_samebases

"""
    lattice_density(state)

Calculate local density of given `state`.
The output of this function is a `LatticeValue`.

## Arguments:
- `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.
"""
function lattice_density(state::StateType{<:OneParticleBasis})
    l = sites(state)
    N = internal_length(state)
    LatticeValue(l, [real(sum(matrix_element(state, j, j) for j in (i-1)*N+1:i*N)) for i in eachindex(l)])
end
function lattice_density(state::StateType{<:ManyBodyBasis{<:OneParticleBasis}})
    vs = zeros(length(basis(state).onebodybasis))
    bas = basis(state)
    for i in 1:length(bas)
        occ = bas.occupations[i]
        vs .+= real.(occ .* matrix_element(state, i, i))
    end
    N = internal_length(state)
    l = sites(state)
    LatticeValue(l, [sum(@view vs[(i-1)*N+1:i*N]) for i in eachindex(l)])
end

function diag_reduce(f, op::OneParticleOperator)
    l = sites(op)
    N = internal_length(op)
    LatticeValue(l,
        [f(@view op.data[(i-1)*N+1:i*N, (i-1)*N+1:i*N]) for i in eachindex(l)])
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
                  for i in eachindex(sites(op)), j in eachindex(sites(op))]
    return AdjacencyMatrix(sites(op), matrix)
end
function adjacency_matrix(op::CompositeLatticeOperator)
    n = internal_length(op)
    ind(k) = (k - 1) * n + 1 : k * n
    matrix = Bool[!iszero(op.data[ind(i), ind(j)])
                  for i in eachindex(sites(op)), j in eachindex(sites(op))]
    return AdjacencyMatrix(sites(op), matrix)
end

function apply_field!(op::AbstractLatticeOperator, field::AbstractField)
    l = sites(op)
    N = internal_length(op)
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            op.data[(i - 1) * N + 1: i * N, (j - 1) * N + 1: j * N] *=
                exp(-2π * im * line_integral(field, site1.coords, site2.coords))
        end
    end
end
