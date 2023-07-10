import QuantumOpticsBase

struct LatticeBasis{LT<:Lattice} <: QuantumOpticsBase.Basis
    shape::SVector{1, Int}
    latt::LT
    LatticeBasis(l::LT) where LT<:Lattice = new{LT}(SA[length(l)], l)
end
Base.:(==)(lb1::LatticeBasis, lb2::LatticeBasis) = lb1.latt == lb2.latt
lattice(lb::LatticeBasis) = lb.latt

QuantumOpticsBase.basisstate(T::Type, b::LatticeBasis, site::LatticeSite) =
    basisstate(T, b, site_index(b.latt, site))
function QuantumOpticsBase.diagonaloperator(lv::LatticeValue)
    diagonaloperator(LatticeBasis(lattice(lv)), lv.values)
end
function QuantumOpticsBase.diagonaloperator(f::Function, b::LatticeBasis)
    diagonaloperator(b, f.(b.latt))
end

function densityoperator(lb::LatticeBasis, l::Lattice)
    check_is_sublattice(lb.latt, l)
    diagonaloperator(in(l), lb)
end

"""
coord_operators(lb::LatticeBasis)

Returns a `Tuple` of coordinate `LatticeOperator`s for given basis.
"""
coord_operators(lb::LatticeBasis) = Tuple(diagonaloperator(lb, lv) for lv in coord_values(lb.latt))
# coord(lb::LatticeBasis, coord) = diagonaloperator(lb, coord_values(lb.latt)[_parse_axis_descriptor(coord)])

const LatticeOperator{T, MT} = Operator{T, T, MT} where T<:LatticeBasis

"""
    apply_field!(hamiltonian, field[; nsteps])

Applies magnetic field to given hamiltonian matrix by adjusting the phase factors.
"""
function apply_field!(ham::LatticeOperator, field::AbstractField)
    l = lattice(ham)
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            if i > j && !iszero(ham[i, j])
                p1 = site1.coords
                p2 = site2.coords
                pmod = exp(-2π * im * line_integral(field, p1, p2))
                !isfinite(pmod) && error("got NaN or Inf when finding the phase factor")
                ham[i, j] *= pmod
                ham[j, i] *= pmod'
            end
        end
    end
end

function add_diagonal!(builder, op, diag)
    for i in 1:length(diag)
        increment!(builder, op * diag[CartesianIndex(i)], i, i)
    end
end

"""
    radius_vector(l::Lattice, hop::Hopping)
Finds the vector between two sites on a lattice according to possibly periodic boundary conditions
(`site2` will be translated along the macrocell to minimize the distance between them).
"""
function radius_vector(l::Lattice, hop::Bonds)
    i, j = hop.site_indices
    bravais(l).basis[:, j] - bravais(l).basis[:, i] +
     mm_assuming_zeros(bravais(l).translation_vectors, hop.translate_uc)
end

@inline _get_bool_value(::Nothing, ::Lattice, ::LatticeSite, ::LatticeSite) = true
@inline _get_bool_value(f::Function, l::Lattice, site1::LatticeSite, site2::LatticeSite) =
    f(l, site1, site2)
@inline _get_bool_value(g::AbstractGraph, ::Lattice, site1::LatticeSite, site2::LatticeSite) =
    match(g, site1, site2)

function add_hoppings!(builder, selector, l::Lattice, op, bond::Bonds,
        field::AbstractField, boundaries::BoundaryConditions)
    dims(bond) > dims(l) && error("Incompatible dims")
    trv = radius_vector(l, bond)
    for site1 in l
        site1.basis_index != hop.site_indices[1] && continue
        site2 = LatticeSite(add_assuming_zeros(site1.unit_cell, hop.translate_uc),
            hop.site_indices[2], site1.coords + trv)
        add_hoppings!(builder, selector, l, op, site1 => site2, field, boundaries)
    end
end

function add_hoppings!(builder, selector, l::Lattice, op, bond::SingleBond,
        field::AbstractField, boundaries::BoundaryConditions)
    site1, site2 = bond
    p1 = site1.coords
    p2 = site2.coords
    factor, site2 = shift_site(boundaries, l, site2)
    i = @inline site_index(l, site1)
    j = @inline site_index(l, site2)
    i === nothing && return
    j === nothing && return
    !_get_bool_value(selector, l, site1, site2) && return
    total_factor = exp(-2π * im * line_integral(field, p1, p2)) * factor
    !isfinite(total_factor) && error("got NaN or Inf when finding the phase factor")
    increment!(builder, op * total_factor, i, j)
    increment!(builder, op' * total_factor', j, i)
end

function hoppings(selector, lb::LatticeBasis, bonds...;
        field::AbstractField=NoField(), boundaries::BoundaryConditions=BoundaryConditions())
    l = lb.latt
    check_lattice_fits(selector, l)
    builder = SparseMatrixBuilder((length(lb), length(lb)))
    for bond in bonds
        add_hoppings!(builder, selector, l, 1, bond, field, boundaries)
    end
    Operator(lb, to_matrix(builder))
end
hoppings(selector, lb::LatticeBasis; kw...) =
    hoppings(selector, lb, default_bonds(lb.latt)...; kw...)
hoppings(selector, l::Lattice, args...; kw...) =
    hoppings(selector, LatticeBasis(l), args...; kw...)

@doc raw"""
    hoppings([f, ]lattice::Lattice, hopping::Hopping[, field::AbstractField])

Creates a hopping operator:
$$\hat{A} = \sum_{pairs} \hat{c}^\dagger_j \hat{c}_i + h. c.$$

Arguments:
- `f`: a function that takes a `Lattice` and two `LatticeSite`s, returns whether this pair should be included.
Can also be a `PairSelector`,
- `lattice`: the lattice to create the operator on.
- `hopping`: the `Hopping` object describing the site pairs and the $\hat{t}$ operator.
- `field`: the `AbstractField` object that defines the magnetic field to generate phase factors using Peierls substitution.
"""
hoppings(l, hops::Bonds...; field::AbstractField=NoField(), boundaries=BoundaryConditions()) =
    hoppings(nothing, l, hops...; field=field, boundaries=boundaries)
