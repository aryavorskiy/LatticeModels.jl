import QuantumOpticsBase: check_samebases

############
# Lattices #
############

"""
    SquareLattice{N}
Type alias for `Lattice{:square,N,1}`.

---
    SquareLattice(sz::Int...)

Constructs a square lattice of size `sz`.
"""
const SquareLattice{N} = Lattice{:square,N,1}
function SquareLattice{N}(sz::Vararg{Int,N}) where {N}
    eye = SMatrix{N,N}(I)
    Lattice(:square, sz, Bravais(eye))
end
default_bonds(::SquareLattice{N}) where {N} = Tuple(SiteOffset(axis=i) for i in 1:N)
default_nnbonds(::SquareLattice{N}) where {N} = Tuple(SiteOffset(one_hot(i, Val(N)) + k * one_hot(j, Val(N))) for i in 1:N for j in 1:i-1 for k in (-1, 1))
default_nnnbonds(::SquareLattice{N}) where {N} = Tuple(SiteOffset(axis=i, dist=2) for i in 1:N)

"""
    HoneycombLattice
Type alias for `Lattice{:honeycomb,2,2}`.

---
    HoneycombLattice(sz::Vararg{Int, 2})

Constructs a honeycomb lattice with a `sz`-size macrocell.
"""
const HoneycombLattice = Lattice{:honeycomb,2,2}
function HoneycombLattice(sz::Vararg{Int, 2})
    bvs = Bravais([1 0.5; 0 √3/2], [0 0.5; 0 √3/6])
    Lattice(:honeycomb, sz, bvs)
end
default_bonds(::HoneycombLattice) = (SiteOffset(2 => 1), SiteOffset(2 => 1, axis=1), SiteOffset(2 => 1, axis=2))
default_nnbonds(::HoneycombLattice) = (
    SiteOffset(1 => 1, axis = 1),
    SiteOffset(2 => 2, axis = 1, dist=-1),
    SiteOffset(1 => 1, axis = 2, dist=-1),
    SiteOffset(2 => 2, axis = 2),
    SiteOffset(1 => 1, [-1, 1]),
    SiteOffset(2 => 2, [1, -1]))
default_nnnbonds(::HoneycombLattice) = (
    SiteOffset(2 => 1, SA[1, 1]),
    SiteOffset(2 => 1, SA[1, -1]),
    SiteOffset(2 => 1, SA[-1, 1]))

##########
# Fields #
##########
"""
    LandauField <: AbstractField

An object representing Landau calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
struct LandauField <: AbstractField
    B::Float64
end
vector_potential(field::LandauField, p1) = (0, p1[1] * field.B)
line_integral(field::LandauField, p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * field.B
Base.show(io::IO, ::MIME"text/plain", field::LandauField) = print(io, "Landau calibration field; B = $(field.B) flux quanta per 1×1 plaquette")

"""
    SymmetricField <: AbstractField

An object representing symmetrically calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
struct SymmetricField <: AbstractField
    B::Float64
end
vector_potential(field::SymmetricField, p1) = SA[-p1[2], p1[1]] * field.B / 2
line_integral(field::SymmetricField, p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * field.B
Base.show(io::IO, ::MIME"text/plain", field::SymmetricField) = print(io, "Symmetric calibration field; B = $(field.B) flux quanta per 1×1 plaquette")

"""
    FluxField <: AbstractField

An object representing a small magnetic flux through given point. The field is directed along z-axis.
Fields:
- `B`: The magnetic field value
- `point`: A `NTuple{2, Number}` representing the point where the magnetic flux is located.
"""
struct FluxField <: AbstractField
    B::Float64
    P::NTuple{2,Float64}
    FluxField(B,P=(0,0)) = new(B, P)
end
function vector_potential(field::FluxField, p1)
    (x, y) = p1
    normsq = (x^2 + y^2)
    SA[-y, x] / normsq * field.B
end
function line_integral(field::FluxField, p1, p2)
    Pv = SVector(field.P)
    p1 = p1[1:2] - Pv
    p2 = p2[1:2] - Pv
    if iszero(p1) || iszero(p2)
        return 0.0
    end
    asin((1 - 1e-11) * det(hcat(p1, p2)) / norm(p1) / norm(p2)) * field.B
end
Base.show(io::IO, ::MIME"text/plain", field::FluxField) = print(io, "Delta flux field through point $(field.P); B = $(field.B) flux quanta")

################
# Hamiltonians #
################

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
qwz(sample::Sample{<:Any, <:SquareLattice}, m=1; kw...) =
    build_hamiltonian(sample,
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
haldane(sample::Sample{<:Any, <:HoneycombLattice}, t1::Real, t2::Real, m::Real=0; kw...) =
    build_hamiltonian(sample,
    (coord_operator(lattice(sample), :index) * 2 - one(LatticeBasis(lattice(sample)))) * m,
    t1 => SiteOffset(2 => 1),
    t1 => SiteOffset(2 => 1, axis = 1),
    t1 => SiteOffset(2 => 1, axis = 2),
    im * t2 => SiteOffset(1 => 1, axis = 1),
    im * t2 => SiteOffset(2 => 2, SA[-1, 0]),
    im * t2 => SiteOffset(1 => 1, SA[0, -1]),
    im * t2 => SiteOffset(2 => 2, axis = 2),
    im * t2 => SiteOffset(1 => 1, SA[-1, 1]),
    im * t2 => SiteOffset(2 => 2, SA[1, -1]); kw...)
@accepts_system haldane

############
# Currents #
############

"""
DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents{HT, DT} <: AbstractCurrents
    hamiltonian::HT
    density::DT

    """
        DensityCurrents(hamiltonian, density_mat)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `density_mat`.
    """
    function DensityCurrents(ham::HT, dens::DT) where {HT<:AbstractLatticeOperator, DT<:AbstractLatticeOperator}
        check_samebases(ham, dens)
        new{HT, DT}(ham, dens)
    end
end

function Base.getindex(curr::DensityCurrents, i::Int, j::Int)
    N = internal_length(curr.hamiltonian)
    is = (i - 1) * N + 1: i * N
    js = (j - 1) * N + 1: j * N
    2imag(tr(curr.density.data[is, js] * curr.hamiltonian.data[js, is]))
end
lattice(curr::DensityCurrents) = lattice(curr.hamiltonian)
