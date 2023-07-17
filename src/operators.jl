import QuantumOpticsBase: basis, samebases, check_samebases

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
coord_operators(l::Lattice) = coord_operators(LatticeBasis(l))
coord(lb::LatticeBasis, coord) = diagonaloperator(lb, [getproperty(site, coord) for site in lb.latt])
coord(l::Lattice, coord) = coord(LatticeBasis(l), coord)

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
hoppings(l, hops::Bonds...; field::AbstractField=NoField(), boundaries=BoundaryConditions()) =
    hoppings(nothing, l, hops...; field=field, boundaries=boundaries)
hoppings(::Nothing, ::Nothing, args...; kw...) =
    throw(MethodError(hoppings, args))

const AbstractBonds = Union{Bonds, SingleBond}
function tight_binding!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractBonds};
        field=NoField(), boundaries=BoundaryConditions())
    # Hopping operator
    opdata, bond = arg
    add_hoppings!(builder, sample.adjacency_matrix, sample.latt, opdata, bond, field, boundaries)
end
function tight_binding!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:LatticeValue}; kw...)
    # Diagonal operator
    opdata, lv = arg
    add_diagonal!(builder, opdata, lv.values)
end
function tight_binding!(builder::SparseMatrixBuilder, ::Sample, arg::Operator; kw...)
    # Arbitrary sparse operator
    increment!(builder, arg.data)
end

function preprocess_argument(sample::Sample, arg::Operator)
    if samebases(basis(arg), onebodybasis(sample))
        sparse(arg)
    elseif samebases(basis(arg), sample.internal)
        sparse(arg) ⊗ one(LatticeBasis(sample.latt))
    elseif samebases(basis(arg), LatticeBasis(sample.latt))
        internal_one(sample) ⊗ sparse(arg)
    else
        error("Invalid Operator argument: basis does not match neither lattice nor internal phase space")
    end
end

function preprocess_argument(sample::Sample, arg::AbstractMatrix)
    bas = sample.internal
    if all(==(length(bas)), size(arg))
        preprocess_argument(sample, SparseOperator(bas, arg))
    else
        error("Invalid Matrix argument: size does not match on-site dimension count")
    end
end

preprocess_argument(sample::Sample, arg) =
    preprocess_argument(sample, internal_one(sample) => arg)

preprocess_argument(sample::Sample, n::Number) = preprocess_argument(sample, n * internal_one(sample))

function preprocess_argument(sample::Sample, arg::Pair)
    op, on_lattice = arg
    if op isa Operator
        check_samebases(basis(op), sample.internal)
        opdata = sparse(op.data)
    elseif op isa AbstractMatrix
        opdata = sparse(op)
    elseif op isa Number
        opdata = op
    else
        error("Invalid Pair argument: unsupported on-site operator type")
    end
    if on_lattice isa LatticeValue
        check_lattice_match(on_lattice, sample)
    elseif !(on_lattice isa AbstractBonds)
        error("Invalid Pair argument: unsupported on-lattice operator type")
    end
    opdata => on_lattice
end

function tight_binding(sample::Sample, args...;
        field=NoField(), boundaries=BoundaryConditions())
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    for arg in args
        tight_binding!(builder, sample, preprocess_argument(sample, arg);
            field=field, boundaries=boundaries)
    end
    bas = onebodybasis(sample)
    oper = Operator(bas, to_matrix(builder))
    if sample.statistics == OneParticle
        return oper
    end
    manybodyoperator(ManyBodyBasis(bas, occupations(sample)), oper)
end
tight_binding(::Nothing, ::Nothing, args...; kw...) = throw(MethodError(tight_binding, args))
