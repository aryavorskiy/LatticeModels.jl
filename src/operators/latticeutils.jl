import QuantumOpticsBase: basis, samebases, check_samebases

"""
    localdensity(state)

Compute the local density of given `state`. The result is a `LatticeValue` with the same
lattice as the input state.

## Arguments
- `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.
"""
function localdensity(state::StateType{<:OneParticleBasis})
    l = lattice(state)
    N = internal_length(state)
    LatticeValue(l, [real(sum(matrix_element(state, j, j) for j in (i-1)*N+1:i*N)) for i in eachindex(l)])
end
function localdensity(state::StateType{<:ManyBodyBasis{<:OneParticleBasis}})
    vs = zeros(length(basis(state).onebodybasis))
    bas = basis(state)
    for i in 1:length(bas)
        occ = bas.occupations[i]
        vs .+= real.(occ .* matrix_element(state, i, i))
    end
    N = internal_length(state)
    l = lattice(state)
    LatticeValue(l, [sum(@view vs[(i-1)*N+1:i*N]) for i in eachindex(l)])
end

function diag_reduce(f, op::OneParticleOperator)
    l = lattice(op)
    N = internal_length(op)
    LatticeValue(l,
        [f(@view op.data[(i-1)*N+1:i*N, (i-1)*N+1:i*N]) for i in eachindex(l)])
end

function QuantumOpticsBase.ptrace(op::CompositeLatticeOperator, sym::Symbol)
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
    adjacencymatrix(op::Operator)

Generates an `AdjacencyMatrix` for the provided operator.
"""
function adjacencymatrix(op::OneParticleOperator)
    n = internal_length(op)
    ind(k) = (k - 1) * n + 1 : k * n
    colptr = Int[]
    rowval = Int[]
    l = lattice(op)
    for i in 1:length(l)
        push!(colptr, length(rowval) + 1)
        for j in 1:length(l)
            if !iszero(op.data[ind(j), ind(i)])
                push!(rowval, j)
            end
        end
    end
    push!(colptr, length(rowval) + 1)
    matrix = SparseMatrixCSC(length(l), length(l), colptr, rowval, fill(true, length(rowval)))
    return AdjacencyMatrix(lattice(op), matrix)
end

function apply_field!(op::OneParticleOperator, field::AbstractField)
    l = lattice(op)
    N = internal_length(op)
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            op.data[(i - 1) * N + 1: i * N, (j - 1) * N + 1: j * N] *=
                exp(-2π * im * line_integral(field, site1.coords, site2.coords))
        end
    end
end
