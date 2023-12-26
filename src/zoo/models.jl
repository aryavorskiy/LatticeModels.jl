hubbard(T::Type, sys::NParticles; U::Real=0, kw...) = tightbinding_hamiltonian(T, sys; kw...) +
    interaction((site1, site2) -> (site1 == site2 ? U : zero(U)), T, sys)
hubbard(t::Type, l::AbstractLattice, N::Int; T = 0, statistics = FermiDirac, kw...) =
    hubbard(t, NParticles(l, N; T = T, statistics = statistics); kw...)
@accepts_t hubbard

@doc raw"""
    qwz([f, ]mv::LatticeValue[; field::AbstractField, pbc=false])
    qwz([f, ]l::SquareLattice[, m::Number=1; field::AbstractField, pbc=false])

$$\hat{H} =
\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i +
\sum_i^\text{sites} \left(
c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i +
c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i +
h. c. \right)$$

Generates a QWZ model hamiltonian operator with set magnetic field.
If the ``m_i`` values are set by the `mv::LatticeValue`, which must be defined on a `SquareLattice`.
Otherwise they will all be set to `m`.

`f` here must be a function or a `PairSelector` describing which hoppings will be excluded.
"""
qwz(m::LatticeValue; kw...) = qwz(lattice(m), m; kw...)
qwz(sys::System{<:Sample{<:SquareLattice}}, m=1; kw...) =
    build_hamiltonian(sys,
    [1 0; 0 -1] => m,
    [1 -im; -im -1] / 2 => SiteOffset(axis = 1),
    [1 -1; 1 -1] / 2 => SiteOffset(axis = 2); kw...)
@accepts_system qwz SpinBasis(1//2)

@doc raw"""
    haldane(l::HoneycombLattice, t1::Real, t2::Real[, m::Real=0; field::AbstractField])

$$\hat{H} =
\sum_i^\text{sublattice A} m c^\dagger_i c_i +
\sum_j^\text{sublattice B} m c^\dagger_j c_j +
\sum_{i, j}^\text{adjacent} \left( t_1 c^\dagger_i c_j + h. c. \right) +
\sum_{i, j}^\text{2-connected,\\counter-clockwise} \left( i \cdot t_2 c^\dagger_i c_j + h. c. \right)$$

Generates a Haldane topological insulator hamiltonian operator.
"""
haldane(sys::System{<:Sample{<:HoneycombLattice}}, t1::Real, t2::Real, m::Real=0; kw...) =
    build_hamiltonian(sys,
    lattice(sys) .|> (site -> site.index == 1 ? m : -m),
    t1 => default_bonds(sys),
    im * t2 => default_bonds(sys, Val(2)); kw...)
@accepts_system haldane

kanemele(sys::System{<:Sample{<:HoneycombLattice}}, t1::Real, t2::Real; kw...) =
    build_hamiltonian(sys,
        t1 => default_bonds(sys),
        im * t2 * sigmaz(internal_basis) => default_bonds(sys, Val(2)); kw...)
@accepts_system kanemele SpinBasis(1//2)
