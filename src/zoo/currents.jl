import QuantumOpticsBase: check_samebases
import QuantumOpticsBase: allocate_buffer, state_transition!, state_index

# Average of creation/destruction operator combination
function _avg(state::StateType{<:OneParticleBasis}, at_indices, a_indices)
    length(at_indices) == 1 && length(a_indices) == 1 || return zero(eltype(state.data))
    at = only(at_indices)
    a = only(a_indices)
    return matrix_element(state, at, a)
end
function _avg(state::StateType{<:ManyBodyBasis{<:OneParticleBasis}}, at_indices, a_indices)
    bas = basis(state)
    buffer = allocate_buffer(bas)
    s = 0.
    for (m, occ) in enumerate(bas.occupations)
        C = state_transition!(buffer, occ, at_indices, a_indices)
        C === nothing && continue
        n = state_index(bas.occupations, buffer)
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
    buffer = allocate_buffer(bas)
    mask = fill(false, N, N)
    for (m, occ) in enumerate(bas.occupations)
        for ni in 1:N, nj in 1:N
            mask[ni, nj] && continue

            i′ = (i - 1) * N + ni
            j′ = (j - 1) * N + nj
            C = state_transition!(buffer, occ, j′, i′)
            C === nothing && continue

            n = state_index(bas.occupations, buffer)
            T[ni, nj] = op.data[m, n] / C
            mask[ni, nj] = true
        end
    end
    return T
end

"""
    DensityCurrents <: AbstractCurrents

Density currents for given state and given hamiltonian.

---
    DensityCurrents(hamiltonian, state)

Constructs a `DensityCurrents` object for given `hamiltonian` and `state`.

## Arguments
- `hamiltonian`: A `Hamiltonian` object representing the Hamiltonian of the system.
- `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.

## Example
## Examples
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(4, 4);

julia> H0 = tightbinding_hamiltonian(lat); psi = groundstate(H0);

julia> H1 = tightbinding_hamiltonian(lat, field=LandauGauge(0.1));

julia> currents = DensityCurrents(H1, psi)
Density currents for system:
One particle on 16-site SquareLattice in 2D space
```
"""
struct DensityCurrents{HT, ST} <: AbstractCurrents
    hamiltonian::HT
    state::ST
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
    for α in 1:N, β in 1:N
        i′ = α + (i - 1) * N
        j′ = β + (j - 1) * N
        outp_cur += 2 * imag(T[α, β] * _avg(curr.state, j′, i′))
    end
    return outp_cur
end

function Base.show(io::IO, mime::MIME"text/plain", curr::DensityCurrents)
    println(io, "Density currents for system:")
    show(io, mime, curr.hamiltonian.sys)
end

"""
    LocalOperatorCurrents <: AbstractCurrents

Local operator (e. g. spin) currents for given state and given hamiltonian.

---
    LocalOperatorCurrents(hamiltonian, state, op)

Constructs a `DensityCurrents` object for given `hamiltonian` and `state`.

## Arguments
- `hamiltonian`: A `Hamiltonian` object representing the Hamiltonian of the system.
- `state`: A `Ket` or `Bra` representing the wavefunction or an `Operator` representing the density matrix.
- `op`: A local (on-site) operator; either an `Operator` or a matrix of such.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(4, 4); site1, site2 = lat[1:2];

julia> H0 = qwz(lat); psi = groundstate(H0);

julia> H1 = qwz(lat, field=LandauGauge(0.1));

julia> op = [1 0; 0 -1];            # Spin operator

julia> currents = LocalOperatorCurrents(H1, psi, op)
Currents of Operator(Spin(1/2))
 1   0
 0  -1
For system:
One particle on (16-site SquareLattice in 2D space) ⊗ Spin(1/2)

julia> op2 = one(SpinBasis(3//2));  # Invalid operator

julia> LocalOperatorCurrents(H1, psi, op2)
ERROR: ArgumentError: Operator must be defined on the internal basis of the Hamiltonian.
[...]
```
"""
struct LocalOperatorCurrents{HT, ST, OT} <: AbstractCurrents
    hamiltonian::HT
    state::ST
    op::OT
    function LocalOperatorCurrents(ham::HT, state::ST, op::OT; check_commutator=true) where {HT<:AbstractLatticeOperator, ST<:StateType, OT<:DataOperator}
        !hasinternal(ham) && throw(ArgumentError("System expected to have internal degrees of freedom"))
        check_samebases(basis(ham), basis(state))
        if internal_basis(ham) == basis(op)
            return new{HT, ST, OT}(ham, state, op)
        else
            throw(ArgumentError("Operator must be defined on the internal basis of the Hamiltonian."))
        end
    end
end
function LocalOperatorCurrents(ham, state, op::AbstractMatrix)
    !hasinternal(ham) && throw(ArgumentError("System expected to have internal degrees of freedom"))
    le = internal_length(ham)
    @check_size op (le, le)
    return LocalOperatorCurrents(ham, state, Operator(internal_basis(ham), op))
end

lattice(curr::LocalOperatorCurrents) = lattice(curr.hamiltonian)
function Base.getindex(curr::LocalOperatorCurrents, i::Int, j::Int)
    outp_cur = 0.
    N = internal_length(curr.hamiltonian)
    T = _block(curr.hamiltonian, i, j)
    O = curr.op.data
    for α in 1:N, β in 1:N
        i′ = α + (i - 1) * N
        j′ = β + (j - 1) * N
        OT_αβ = sum(O[α, k] * T[k, β] for k in 1:N)
        outp_cur += 2 * imag(OT_αβ * _avg(curr.state, j′, i′))
    end
    return outp_cur
end

function Base.show(io::IO, mime::MIME"text/plain", curr::LocalOperatorCurrents)
    print(io, "Currents of Operator(")
    show(io, mime, basis(curr.op))
    println(io, ")")
    Base.print_array(io, curr.op.data)
    println(io, "\nFor system:")
    show(io, mime, curr.hamiltonian.sys)
end
