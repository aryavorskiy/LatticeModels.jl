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
default_bonds(::SquareLattice{N}) where {N} = Tuple(Bonds(axis=i) for i in 1:N)

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
default_bonds(::HoneycombLattice) = (Bonds(2 => 1), Bonds(2 => 1, axis=1), Bonds(2 => 1, axis=2))

##########
# Fields #
##########

@field_def struct LandauField(B::Number)
    vector_potential(x) = (0, x*B)
    line_integral(p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * B
    show(io::IO, ::MIME"text/plain") = print(io, "Landau calibration field; B = $B flux quanta per 1×1 plaquette")
end
"""
    LandauField <: AbstractField

An object representing Landau calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
LandauField

@field_def struct SymmetricField(B::Number)
    vector_potential(x, y) = SA[-y, x] * B / 2
    line_integral(p1, p2) = (p1[1] * p2[2] - p2[1] * p1[2]) / 2 * B
    show(io::IO, ::MIME"text/plain") = print(io, "Symmetric calibration field; B = $B flux quanta per 1×1 plaquette")
end
"""
    SymmetricField <: AbstractField

An object representing symmetrically calibrated uniform magnetic field along z-axis.
Fields:
- `B`: The magnetic field value
"""
SymmetricField

_angle(p1, p2) = asin((1 - 1e-11) * det(hcat(p1, p2)) / norm(p1) / norm(p2))
@field_def struct FluxField(B::Number, P::NTuple{2,Number} = (0, 0))
    function vector_potential(x, y)
        norm = (x^2 + y^2)
        (-y / norm * B, x / norm * B)
    end
    function line_integral(p1, p2)
        Pv = SVector(P)
        p1 = p1[1:2] - Pv
        p2 = p2[1:2] - Pv
        if iszero(p1) || iszero(p2)
            return 0.0
        end
        _angle(p1, p2) * B
    end
    show(io::IO, ::MIME"text/plain") = print(io, "Delta flux field through point $P; B = $B flux quanta")
end
"""
    FluxField <: AbstractField

An object representing a small magnetic flux through given point. The field is directed along z-axis.
Fields:
- `B`: The magnetic field value
- `point`: A `NTuple{2, Number}` representing the point where the magnetic flux is located.
"""
FluxField

################
# Hamiltonians #
################

# TODO: fix these

# @doc raw"""
#     qwz([f, ]mv::LatticeValue[; field::AbstractField, pbc=false])
#     qwz([f, ]l::SquareLattice[, m::Number=1; field::AbstractField, pbc=false])

# $$\hat{H} =
# \sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i +
# \sum_i^\text{sites} \left(
# c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i +
# c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i +
# h. c. \right)$$

# Generates a spin-orbital tight-binding hamiltonian operator with set magnetic field and boundary conditions.
# If the ``m_i`` values are set by the `mv::LatticeValue`, which must be defined on a `SquareLattice`.
# Otherwise they will all be set to `m`.

# `f` here must be a function or a `PairSelector` describing which hoppings will be excluded.
# """
# SpinTightBinding(f, m::LatticeValue{<:Number, :square}; field=NoField(), pbc=false) =
# @hamiltonian begin
#     lattice := lattice(m)
#     dims_internal := 2
#     field := field
#     @diag m ⊗ [1 0; 0 -1]
#     @hop axis=1 [1 -im; -im -1] / 2 pbc=pbc f
#     @hop axis=2 [1 -1; 1 -1] / 2 pbc=pbc f
# end
# SpinTightBinding(f, l::SquareLattice, m::Number=1; field=NoField(), pbc=false) =
# @hamiltonian begin
#     lattice := l
#     dims_internal := 2
#     field := field
#     @diag [m 0; 0 -m]
#     @hop axis=1 [1 -im; -im -1] / 2 pbc=pbc f
#     @hop axis=2 [1 -1; 1 -1] / 2 pbc=pbc f
# end
# SpinTightBinding(f, l_sz::NTuple{N, Int}, m::Number=1; kw...) where N =
#     SpinTightBinding(f, SquareLattice(l_sz...), m; kw...)
# SpinTightBinding(args...; kw...) = SpinTightBinding(nothing, args...; kw...)
# SpinTightBinding(::Nothing, ::Nothing, args...; kw...) = throw(MethodError(SpinTightBinding, args))

@doc raw"""
    haldane(l::HoneycombLattice, t1::Real, t2::Real[, m::Real=0; field::AbstractField])

$$\hat{H} =
\sum_i^\text{sublattice A} m c^\dagger_i c_i +
\sum_j^\text{sublattice B} m c^\dagger_j c_j +
\sum_{i, j}^\text{adjacent} \left( t_1 c^\dagger_i c_j + h. c. \right) +
\sum_{i, j}^\text{2-connected,\\counter-clockwise} \left( i \cdot t_2 c^\dagger_i c_j + h. c. \right)$$

Generates a Haldane topological insulator hamiltonian operator.
"""
haldane(l::HoneycombLattice, t1::Real, t2::Real, m::Real=0; kw...) =
    tight_binding(l,
    (coord(l, :basis_index) * 2 - one(LatticeBasis(l))) * m,
    t1 => Bonds(2 => 1),
    t1 => Bonds(2 => 1, axis = 1),
    t1 => Bonds(2 => 1, axis = 2),
    im * t2 => Bonds(1 => 1, axis = 1),
    -im * t2 => Bonds(2 => 2, axis = 1),
    -im * t2 => Bonds(1 => 1, axis = 2),
    im * t2 => Bonds(2 => 2, axis = 2),
    im * t2 => Bonds(1 => 1, [-1, 1]),
    -im * t2 => Bonds(2 => 2, [-1, 1]);
    kw...)

############
# Currents #
############

"""
DensityCurrents <: AbstractCurrents

Density currents for given density matrix and given hamiltonian.
"""
struct DensityCurrents <: AbstractCurrents
    hamiltonian::LatticeOperator
    density::LatticeOperator

    """
        DensityCurrents(hamiltonian, density_mat)

    Constructs a `DensityCurrents` object for given `hamiltonian` and `density_mat`.
    """
    function DensityCurrents(ham::LatticeOperator, dens::LatticeOperator)
        check_basis_match(ham, dens)
        new(ham, dens)
    end
end

Base.getindex(curr::DensityCurrents, i::Int, j::Int) =
    2imag(tr(curr.density[i, j] * curr.hamiltonian[j, i]))
lattice(curr::DensityCurrents) = curr.hamiltonian.basis.lattice

==#
