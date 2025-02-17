"""
    hubbard([type, ]sys[; U, T, t1, t2, t3, field])

Generates a Hubbard model hamiltonian on given manybody system `sys`.

See [`fermihubbard`](@ref), [`bosehubbard`](@ref) for more specific models.
"""
function hubbard(T::Type, sys::ManyBodySystem, args...; U::Real=0, kw...)
    if sys isa NParticles && sys.statistics == FermiDirac && !hasinternal(sys)
        @warn """Fermi-Dirac statistics with no internal degrees of freedom.
        No interaction is possible."""
    end
    op = Operator(tightbinding_hamiltonian(T, sys, args...; kw...)) +
        interaction((site1, site2) -> (site1 == site2 ? U : zero(U)), T, sys)
    return Hamiltonian(sys, op)
end
"""
    bosehubbard([type, ]lat, N[; U, T, t1, t2, t3, field])

``\\hat{H} =
\\sum_{i,j}^\\text{sites} t_{ij} c^\\dagger_i c_j +
\\sum_i^\\text{sites} \\frac{U}{2} \\hat{n}_i (\\hat{n}_i - 1)``

Generates a Bose-Hubbard model hamiltonian on given lattice `lat`.
## Arguments
- `type`: The element type of the resulting operator. Default is `ComplexF64`.
- `N`: The number of particles.

## Keyword arguments
- `t1`, `t2`, `t3` denote the coefficient on first, second and third hoppings respectively.
    By default `t1` is equal to one, the rest are zero.
- `U`: The interaction strength. Default is zero.
- `T`: The temperature of the system. Default is zero.
- `field`: The magnetic field. Default is `NoField()`.
"""
bosehubbard(type::Type, l::AbstractLattice, N::Int; T = 0, occupations_type=nothing, kw...) =
    hubbard(type, NParticles(l, N; T = T, statistics = BoseEinstein, occupations_type=occupations_type); kw...)

_nspins(v::AbstractVector, up=true) = sum(@view v[2 - up:2:end])
_odds_mask(v::Unsigned) = typemax(v) ÷ 3
function _odds_mask(v::BigInt)
    lenp = prevpow(2, v)
    mask = lenp + (lenp - 1) ÷ 3
    if mask & 1 == 0
        mask >>= 1
    end
    return mask
end
function _nspins(v::FermionBitstring, up=true)
    mask = _odds_mask(v.bits) << (1 - up)
    return count_ones(v.bits & mask)
end
function _filter_occupations!(occ::AbstractVector, N_up, N_down)
    keep_inds = Int[]
    for i in eachindex(occ)
        if _nspins(occ[i]) != N_up || _nspins(occ[i], false) != N_down
            push!(keep_inds, i)
        end
    end
    if occ isa QuantumOpticsBase.SortedVector
        deleteat!(occ.sortedvector, keep_inds)
    else
        deleteat!(occ, keep_inds)
    end
end

"""
    FermiHubbardSpinSystem(l, N_up, N_down[; occupations_type])

Generates a Fermi-Hubbard model basis system on given lattice `l`. The number of particles is
`N_up` and `N_down` for spin up and down respectively. The `occupations_type` keyword argument
defines the type of the occupation basis.
"""
function FermiHubbardSpinSystem(l::AbstractLattice, N_up, N_down; occupations_type=nothing, kw...)
    occ = occupations(NParticles(l ⊗ SpinBasis(1//2), N_up + N_down, occupations_type=occupations_type))
    _filter_occupations!(occ, N_up, N_down)
    bas = ManyBodyBasis(basis(l ⊗ SpinBasis(1//2)), occ)
    return ManyBodyBasisSystem(bas; kw...)
end

"""
    fermihubbard([type, ]lat, N[; U, T, t1, t2, t3, field])

``\\hat{H} =
\\sum_{i,j}^\\text{sites} t_{ij} c^\\dagger_i c_j +
\\sum_i^\\text{sites} U \\hat{n}_i^{\\uparrow} \\hat{n}_i^{\\downarrow}``

Generates a Fermi-Hubbard model hamiltonian on given lattice `lat`.

## Arguments
- `type`: The element type of the resulting operator. Default is `ComplexF64`.
- `N`: The number of particles.

## Keyword arguments
- `t1`, `t2`, `t3` denote the coefficient on first, second and third hoppings respectively.
    By default `t1` is equal to one, the rest are zero.
- `U`: The interaction strength. Default is zero.
- `T`: The temperature of the system. Default is zero.
- `field`: The magnetic field. Default is `NoField()`.
"""
fermihubbard(type::Type, l::AbstractLattice, N::Int; T = 0, occupations_type=nothing, kw...) =
    hubbard(type, NParticles(l ⊗ SpinBasis(1//2), N; T = T, statistics = FermiDirac, occupations_type=occupations_type); kw...)
@accepts_t hubbard
@accepts_t bosehubbard
@accepts_t fermihubbard

"""
    qwz(m[; T, μ, field, statistics])
    qwz(lat[, m; T, μ, field, statistics])

``\\hat{H} =
\\sum_i^\\text{sites} m_i c^\\dagger_i \\sigma_z c_i +
\\sum_i^\\text{sites} \\left(
c^\\dagger_{i + \\hat{x}} \\frac{\\sigma_z - i \\sigma_x}{2} c_i +
c^\\dagger_{i + \\hat{y}} \\frac{\\sigma_z - i \\sigma_y}{2} c_i +
h. c. \\right)``

Generates a QWZ model hamiltonian operator on given square lattice `lat`.

## Arguments
- `m` (either a `LatticeValue` or a number) defines the ``m_i`` factors

## Keyword arguments
- `T`: The temperature of the system. Default is zero.
- `μ`: The chemical potential. Use `mu` as a synonym if Unicode input is not available.
- `field`: The magnetic field. Default is `NoField()`.
- `statistics` defines the particle statistics, either `FermiDirac` or `BoseEinstein`.
"""
qwz(m::LatticeValue; kw...) = qwz(lattice(m), m; kw...)
function qwz(sys::System, m=1; kw...)
    checklatticetype(sys, SquareLattice)
    construct_hamiltonian(sys,
    [1 0; 0 -1] => m,
    [1 -im; -im -1] / 2 => BravaisTranslation(axis = 1),
    [1 -1; 1 -1] / 2 => BravaisTranslation(axis = 2); kw...)
end
@accepts_system qwz SpinBasis(1//2)

const honeycomb_2nn = BravaisSiteMapping(
    BravaisTranslation(1 => 1, axis = 1),
    BravaisTranslation(2 => 2, axis = 1, dist=-1),
    BravaisTranslation(1 => 1, axis = 2, dist=-1),
    BravaisTranslation(2 => 2, axis = 2),
    BravaisTranslation(1 => 1, [-1, 1]),
    BravaisTranslation(2 => 2, [1, -1]))

"""
    haldane(lat, t1, t2[, m=0; T, μ, field, statistics])

``\\hat{H} =
\\sum_i^\\text{sublattice A} m c^\\dagger_i c_i +
\\sum_j^\\text{sublattice B} m c^\\dagger_j c_j +
\\sum_{i, j}^\\text{adjacent} \\left( t_1 c^\\dagger_i c_j + h. c. \\right) +
\\sum_{i, j}^\\text{2-connected, counter-clockwise} \\left( i \\cdot t_2 c^\\dagger_i c_j + h. c. \\right)``

Generates a Haldane topological insulator hamiltonian operator on given lattice `lat`.

## Keyword arguments
- `T`: The temperature of the system. Default is zero.
- `μ`: The chemical potential. Use `mu` as a synonym if Unicode input is not available.
- `field`: The magnetic field. Default is `NoField()`.
- `statistics` defines the particle statistics, either `FermiDirac` or `BoseEinstein`.
"""
function haldane(sys::System, t1::Real, t2::Real, m::Real=0; kw...)
    checklatticetype(sys, HoneycombLattice)
    construct_hamiltonian(sys,
    lattice(sys) .|> (site -> site.basindex == 1 ? m : -m),
    t1 => NearestNeighbor(1),
    im * t2 => honeycomb_2nn; kw...)
end
@accepts_system haldane

"""
    kanemele(lat, t1, t2[; T, μ, field, statistics])

``\\hat{H} =
\\sum_{i, j}^\\text{adjacent} \\left( t_1 c^\\dagger_i c_j + h. c. \\right) +
\\sum_{i, j}^\\text{2-connected, counter-clockwise} \\left( i \\cdot t_2 c^\\dagger_i σ_z c_j + h. c. \\right)``

Generates a Kane-Mele hamiltonian operator on given lattice `lat`.

## Keyword arguments
- `T`: The temperature of the system. Default is zero.
- `μ`: The chemical potential. Use keyword `mu` as a synonym if Unicode input is not available.
- `field`: The magnetic field. Default is `NoField()`.
- `statistics` defines the particle statistics, either `FermiDirac` or `BoseEinstein`.
"""
function kanemele(sys::System, t1::Real, t2::Real; kw...)
    checklatticetype(sys, HoneycombLattice)
    construct_hamiltonian(sys,
        t1 => NearestNeighbor(1),
        im * t2 * sigmaz(internal_basis(sys)) => honeycomb_2nn; kw...)
end
@accepts_system kanemele SpinBasis(1//2)
