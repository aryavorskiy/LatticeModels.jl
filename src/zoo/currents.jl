import QuantumOpticsBase: check_samebases

# Average of creation/destruction operator combination
function _avg(state::StateType{<:OneParticleBasis}, at_indices, a_indices)
    length(at_indices) == 1 && length(a_indices) == 1 || return zero(eltype(state.data))
    at = only(at_indices)
    a = only(a_indices)
    return matrix_element(state, a, at)
end
function _avg(state::StateType{<:ManyBodyBasis{<:OneParticleBasis}}, at_indices, a_indices)
    bas = basis(state)
    buffer = QuantumOpticsBase.allocate_buffer(bas.occupations)
    s = 0.
    for (m, occ) in enumerate(bas.occupations)
        C = QuantumOpticsBase.state_transition!(buffer, occ, at_indices, a_indices)
        C === nothing && continue
        n = QuantumOpticsBase.state_index(bas.occupations, buffer)
        s += C * matrix_element(state, n, m)
    end
    return s
end

function _block(op::OneParticleOperator, i, j)
    N = internal_length(op)
    return @view op.data[(i - 1) * N + 1:i * N, (j - 1) * N + 1:j * N]
end
function _block(op::AbstractLatticeOperator, i, j)
    N = internal_length(op)
    T = zeros(eltype(op.data), N, N)
    bas = basis(op)
    buffer = QuantumOpticsBase.allocate_buffer(bas.occupations)
    mask = fill(false, N, N)
    for (m, occ) in enumerate(bas.occupations)
        for ni in 1:N, nj in 1:N
            mask[ni, nj] && continue

            i′ = (i - 1) * N + ni
            j′ = (j - 1) * N + nj
            C = QuantumOpticsBase.state_transition!(buffer, occ, i′, j′)
            C === nothing && continue

            n = QuantumOpticsBase.state_index(bas.occupations, buffer)
            T[ni, nj] = op.data[m, n] / C
            mask[ni, nj] = true
        end
    end
    return T
end

"""
DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents{HT, ST} <: AbstractCurrents
    hamiltonian::HT
    state::ST

    """
        DensityCurrents(hamiltonian, state)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `state`.

    ## Arguments:
    - `hamiltonian`: A `Hamiltonian` object representing the Hamiltonian of the system.
    - `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.
    """
    function DensityCurrents(ham::HT, state::ST) where {HT<:AbstractLatticeOperator, ST<:StateType}
        check_samebases(basis(ham), basis(state))
        new{HT, ST}(ham, state)
    end
end

lattice(curr::DensityCurrents) = lattice(curr.hamiltonian)

function Base.getindex(curr::DensityCurrents, i::Int, j::Int)
    outp_cur = 0.
    N = internal_length(curr.hamiltonian)
    T = _block(curr.hamiltonian, i, j)
    for i′ in 1:N, j′ in 1:N
        i′′ = i′ + (i - 1) * N
        j′′ = j′ + (j - 1) * N
        outp_cur -= 2 * imag(T[i′, j′] * _avg(curr.state, i′′, j′′))
    end
    return outp_cur
end

struct OperatorCurrents{HT, ST, OT} <: AbstractCurrents
    hamiltonian::HT
    state::ST
    op::OT

    """
    LocalOperatorCurrents(hamiltonian, state, op)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `state`.

    ## Arguments:
    - `hamiltonian`: A `Hamiltonian` object representing the Hamiltonian of the system.
    - `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.
    - `op`: A local (on-site) operator; either an `Operator` or a matrix of such.
    """
    function OperatorCurrents(ham::HT, state::ST, op::OT; check_commutator=true) where {HT<:AbstractLatticeOperator, ST<:StateType, OT<:DataOperator}
        !hasinternal(ham) && throw(ArgumentError("System expected to have internal degrees of freedom"))
        check_samebases(basis(ham), basis(state))
        if basis(ham) == basis(op)
            if check_commutator
                comm = ham * op - op * ham
                isapprox(comm.data, zero(comm).data, atol=1e-10) || @warn "The operator must commute with the Hamiltonian"
            end
            return new{HT, ST, OT}(ham, state, op)
        elseif hasinternal(ham) && internal_basis(ham) == basis(op)
            return new{HT, ST, OT}(ham, state, op)
        end
    end
end
function OperatorCurrents(ham, state, op::AbstractMatrix)
    !hasinternal(ham) && throw(ArgumentError("System expected to have internal degrees of freedom"))
    new_op = Operator(internal_basis(ham), op)
    return OperatorCurrents(ham, state, new_op)
end

_op_diag_block(op::DataOperator, _) = op.data
_op_diag_block(op::LatticeOperator, i) = _block(op, i, i)
lattice(curr::OperatorCurrents) = lattice(curr.hamiltonian)
function Base.getindex(curr::OperatorCurrents, i::Int, j::Int)
    outp_cur = 0.
    N = internal_length(curr.hamiltonian)
    T = _block(curr.hamiltonian, i, j)
    O = _op_diag_block(curr.op, i)
    for i′ in 1:N, j′ in 1:N
        i′′ = i′ + (i - 1) * N
        j′′ = j′ + (j - 1) * N
        OT_ij = sum(O[i′, k] * T[k, j′] for k in 1:N)
        TO_ij = sum(T[i′, k] * O[k, j′] for k in 1:N)
        outp_cur -= 2 * imag(_avg(curr.state, i′′, j′′) * OT_ij)
        # TODO: what the hell?
        # outp_cur -= 2 * imag(_avg(curr.state, (i′′, i′′), (i′′, j′′)) * (TO_ij - OT_ij))
    end
    return outp_cur
end
