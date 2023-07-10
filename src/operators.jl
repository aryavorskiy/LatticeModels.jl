include("operators_core.jl")

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
coord(lb::LatticeBasis, coord) = diagonaloperator(lb, [getproperty(site, coord) for site in lb.latt])

function hoppings(selector, lb::LatticeBasis, bonds...;
        field::AbstractField=NoField(), boundaries::BoundaryConditions=BoundaryConditions())
    l = lb.latt
    check_lattice_fits(selector, l)
    builder = SparseMatrixBuilder{ComplexF64}(length(lb), length(lb))
    for bond in bonds
        add_hoppings!(builder, selector, l, 1, bond, field, boundaries)
    end
    Operator(lb, to_matrix(builder))
end
hoppings(selector, lb::LatticeBasis; kw...) =
    hoppings(selector, lb, default_bonds(lb.latt)...; kw...)
hoppings(selector, l::Lattice, args...; kw...) =
    hoppings(selector, LatticeBasis(l), args...; kw...)

const AbstractBonds = Union{Bonds, SingleBond}
function tight_binding!(selector, builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractBonds};
        field=NoField(), boundaries=BoundaryConditions())
    op, bond = arg
    opdata = op isa Operator ? op.data : op
    add_hoppings!(builder, selector, sample.latt, opdata, bond, field, boundaries)
end

function tight_binding!(selector, builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:LatticeValue}; kw...)
    op, lv = arg
    check_lattice_match(sample, lv)
    opdata = op isa Operator ? op.data : op
    add_diagonal!(builder, opdata, lv.values)
end

function process_argument(sample, arg)
    internal_one(sample) => arg
end

function process_argument(sample::Sample, arg::Operator)
    if QuantumOpticsBase.samebases(arg, sample.internal)
        arg.data => ones(sample.latt)
    elseif QuantumOpticsBase.samebases(arg, LatticeBasis(sample.latt))
        internal_one(sample).data => arg.data
    else
        error("Invalid Operator argument: basis does not match neither lattice nor internal phase space")
    end
end

function process_argument(sample::Sample, arg::Pair)
    op, on_lattice = arg
    if op isa Operator
        QuantumOpticsBase.check_samebases(QuantumOpticsBase.basis(op), sample.internal)
        opdata = op.data
    else
        opdata = op
    end
    opdata => on_lattice
end

function tight_binding(selector, sample::Sample, args...;
        field=NoField(), boundaries=BoundaryConditions())
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    for arg in args
        tight_binding!(selector, builder, sample, process_argument(sample, arg);
            field=field, boundaries=boundaries)
    end
    bas = one_particle_basis(sample)
    oper = Operator(bas, to_matrix(builder))
    if sample.statistics == one_particle
        return oper
    end
    occupations = sample.statistics == fermi ? fermionstates(bas, sample.nparticles) : bosonstates(bas, sample.nparticles)
    manybodyoperator(ManyBodyBasis(bas, occupations), oper)
end

tight_binding(selector, latt::Lattice, args...; kw...) = tight_binding(selector, Sample(latt), args...; kw...)
tight_binding(sample, args...; kw...) = tight_binding(nothing, sample, args...; kw...)
tight_binding(::Nothing, ::Nothing, args...; kw...) = throw(MethodError(tight_binding, args))

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
