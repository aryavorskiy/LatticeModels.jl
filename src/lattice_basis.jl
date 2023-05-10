"""
    SingleParticleBasis{LT} where {LT<:Lattice}

A basis on a lattice with some number of internal states on each site.
Fields:
- `lattice`: the [`Lattice`](@ref) of the basis
- `internal_dim`: the number of internal states on each site
"""
struct SingleParticleBasis{LT<:Lattice} <: Basis
    lattice::LT
    internal_dim::Int
end
Basis(l::Lattice, N::Int) = SingleParticleBasis(l, N)
lattice(b::SingleParticleBasis) = b.lattice
dims_internal(b::SingleParticleBasis) = b.internal_dim
==(b1::SingleParticleBasis, b2::SingleParticleBasis) =
    b1.internal_dim == b2.internal_dim && b1.lattice == b2.lattice
site_states(b::SingleParticleBasis) = lattice(b)

length(b::SingleParticleBasis) = length(lattice(b)) * dims_internal(b)
to_slice(b::SingleParticleBasis, i::Int) = b.internal_dim * (i - 1) + 1:b.internal_dim * i
to_slice(b::SingleParticleBasis, site::LatticeSite) = to_slice(b, site_index(lattice(b), site))

function show(io::IO, m::MIME"text/plain", b::SingleParticleBasis)
    println(io, "Single-particle basis with $(b.internal_dim) internal degrees of freedom")
    print(io, "on ")
    show(io, m, b.lattice)
end

const STATS_VARIANTS = (:fermi, :bose)
struct MultiParticleBasis{ST, LT<:Lattice} <: Basis
    lattice::LT
    n_particles::Int
    function MultiParticleBasis(l::Lattice; N, st)
        st ∉ STATS_VARIANTS && error("Insupported statistics :$st, expected one of $STATS_VARIANTS")
        new{st, typeof(l)}(l, N)
    end
end
Basis(l::Lattice; N::Int, st::Symbol=:bose) = MultiParticleBasis(l, N=N, st=st)
lattice(b::MultiParticleBasis) = b.lattice
dims_internal(b::MultiParticleBasis) = 1
==(b1::T, b2::T) where T<:MultiParticleBasis =
    b1.n_particles == b2.particle && b1.lattice == b2.lattice
site_states(b::MultiParticleBasis) = b

length(b::MultiParticleBasis{:bose}) = binomial(length(b.lattice) + b.n_particles - 1, b.n_particles)
length(b::MultiParticleBasis{:fermi}) = binomial(length(b.lattice), b.n_particles)
to_slice(::MultiParticleBasis, i::Int) = i:i

function show(io::IO, m::MIME"text/plain", b::MultiParticleBasis{ST}) where ST
    println(io, "$(b.n_particles)-particle basis with $ST statistics")
    print(io, "on ")
    show(io, m, b.lattice)
end

struct MultiParticleState{LT}
    l::LT
    fills::IdDict{Int,Int}
end
getindex(st::MultiParticleState, i) = get(st.fills, i, 0)
function setindex!(st::MultiParticleState, rhs, i)
    rhs == 0 ? pop!(st.fills, i) : setindex!(st.fills, rhs, i)
    return rhs
end
pairs(st::MultiParticleState) = Iterators.map(st.fills) do (i, n)
    st.l[i] => n
end
length(m::MultiParticleState) = sum(values(m.fills))
function iterate(m::MultiParticleState)
    ks = keys(m.fills)
    istate = iterate(ks)
    istate === nothing && return nothing
    index, iter = istate
    site = m.l[index]
    site, (site, 2, m.fills[index], index, iter)
end
function iterate(m::MultiParticleState, (site, r, f, index, iter))
    r ≤ f && return site, (site, r+1, f, index, iter)
    istate = iterate(keys(m.fills), iter)
    istate === nothing && return nothing
    index, iter = istate
    site = m.l[index]
    site, (site, 2, m.fills[index], index, iter)
end

function iterate(b::MultiParticleBasis{:bose})
    st = MultiParticleState(b.lattice, IdDict(1=>b.n_particles))
    (st, (st, 1, length(b.lattice)))
end
function iterate(::MultiParticleBasis{:bose}, (st, bottom, len))
    bottom == len && return nothing
    va = st[bottom]
    if bottom == 1
        st[1] = va - 1
        st[2] += 1
        va == 1 && (bottom = 2)
    elseif va == 1
        st[bottom] = 0
        st[bottom + 1] += 1
        bottom += 1
    else
        st[bottom] = 0
        st[bottom + 1] += 1
        st[1] = va - 1
        bottom = 1
    end
    st, (st, bottom, len)
end

"""
    TensorProduct{LVT, MT} where {LVT<:LatticeValue{<:Number}, MT<:AbstractMatrix}

A lazy representation of an operator as a tensor product of two distinct phase spaces.
One affects only the internal space, the other - only the lattice space.

The `lattice_value ⊗ matrix` notation computes the value of the `TensorProduct` eagerly,
which means that the result will be a `LatticeOperator`.
However, in the `@hamiltonian` macro lazy computation is forced.
"""
struct TensorProduct{LVT<:LatticeValue{<:Number},N,T}
    lattice_value::LVT
    matrix::SMatrix{N,N,T}
    function TensorProduct(lv::LVT, m::AbstractMatrix{T}) where {LVT,T}
        N = size(m)[1]
        new{LVT,N,T}(lv, m)
    end
end

dims_internal(tp::TensorProduct) = size(tp.matrix)[1]
lattice(tp::TensorProduct) = lattice(tp.lattice_value)
basis(tp::TensorProduct) = Basis(lattice(tp), dims_internal(tp))
zero(tp::TensorProduct) = zero_on_basis(lattice(tp), tp.matrix)
materialize(tp::TensorProduct) = _diag_operator!(zero_on_basis(basis(tp)), tp)
⊗(lv::LatticeValue, m::Matrix) = materialize(TensorProduct(lv, m))
⊗(m::Matrix, lv::LatticeValue) = materialize(TensorProduct(lv, m))

zero_on_basis(bas::Basis) = LatticeArray(bas, zeros(ComplexF64, length(bas), length(bas)))
function zero_on_basis(bas::Basis, MT::Type{<:AbstractMatrix})
    arr = similar(MT, (length(bas), length(bas)))
    LatticeArray(bas, fill!(arr, 0))
end
zero_on_basis(bas::Basis, ::Type{SparseMatrixBuilder{T}}) where T =
    LatticeArray(bas, SparseMatrixBuilder{T}((length(bas), length(bas))))

_wrap_eye(n::Number, eye::Matrix) = n * eye
_wrap_eye(m::AbstractMatrix, ::Matrix) = m
_wrap_eye(::T, ::Matrix) where T = error("Lambda returned a $T, expected Number or AbstractMatrix")
@inline _get_matrix_value(f::Function, state, ::Int, eye::Matrix) = _wrap_eye(f(state), eye)
@inline _get_matrix_value(m::AbstractMatrix, _, ::Int, ::Matrix) = m
@inline _get_matrix_value(tp::TensorProduct, _, i::Int, ::Matrix) = tp.lattice_value.values[i] * tp.matrix
@inline _get_matrix_value(lv::LatticeValue, _, i::Int, eye::Matrix) = lv.values[i] * eye
@inline _get_matrix_value(n::Number, _, i::Int, eye::Matrix) = n * eye
function _diag_operator!(lop::LatticeOperator, op_object)
    N = dims_internal(lop)
    eye = Matrix(I, N, N)
    for (i, state) in enumerate(site_states(lop))
        increment!(lop, _get_matrix_value(op_object, state, i, eye), i, i)
    end
    lop
end

"""
    diag_operator(f, bas::Basis)
    diag_operator(f, l::Lattice, N::Int)

Creates a diagonal operator by applying the `f` function to each site of the lattice of given basis.
`f` must accept a `LatticeSite` and its coordinate vector and return a number or a matrix
which represents operator affecting the internal state of the site.
"""
diag_operator(f::Function, bas::Basis) = _diag_operator!(zero_on_basis(bas), f)

"""
    diag_operator(lv::LatticeValue, N::Int=1)

Creates a diagonal operator which affects only the lattice space.
The `lv` argument must be a `LatticeValue` storing diagonal elements of the operator in lattice space.
`N` is the number of internal degrees of freedom on each site.
"""
diag_operator(lv::LatticeValue{<:Number}, N::Int=1) = lv ⊗ Matrix(I, N, N)

"""
    coord_operators(basis::Basis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
function coord_operators(bas::SingleParticleBasis)
    N = dims_internal(bas)
    d = dims(bas)
    eye = Matrix(I, N, N)
    xyz_operators = [LatticeArray(bas, op_mat) for op_mat in
                     eachslice(zeros(length(bas), length(bas), d), dims=3)]
    for (i, site) in enumerate(lattice(bas))
        for j in 1:d
            xyz_operators[j][i, i] = site.coords[j] * eye
        end
    end
    xyz_operators
end

"""
    coord_operators(lattice::Lattice, ndims::Tnt)

The same as `coord_operators(Basis(lattice, ndims))`.
"""
coord_operators(l::Lattice, N::Int) = coord_operators(Basis(l, N))
