import QuantumOpticsBase: DataOperator

function add_diagonal!(builder, op, diag)
    for i in eachindex(diag)
        increment!(builder, op, i, i, factor=diag[CartesianIndex(i)])
    end
end

function add_hoppings!(builder, l::AbstractLattice, op, bond::SiteOffset,
        field::AbstractField)
    dims(bond) > dims(l) && error("Incompatible dims")
    for site1 in l
        site2 = site1 + bond
        site2 === nothing && continue
        add_hoppings!(builder, l, op, site1 => site2, field)
    end
end

function add_hoppings!(builder, l::AbstractLattice, op, bond::SingleBond,
        field::AbstractField)
    site1, site2 = bond
    p1 = site1.coords
    p2 = site2.coords
    factor, site2 = shift_site(l, site2)
    i = site_index(l, site1)
    j = site_index(l, site2)
    i === nothing && return
    j === nothing && return
    total_factor = exp(-2π * im * line_integral(field, p1, p2)) * factor
    !isfinite(total_factor) && error("got NaN or Inf when finding the phase factor")
    increment!(builder, op, i, j, factor=total_factor)
    increment!(builder, op', j, i, factor=total_factor')
end

struct Hamiltonian{SystemT, BasisT, T} <: DataOperator{BasisT, BasisT}
    sys::SystemT
    basis_l::BasisT
    basis_r::BasisT
    data::T
end
function Hamiltonian(sys::System, op::Operator)
    return Hamiltonian(sys, basis(op), basis(op), op.data)
end
QuantumOpticsBase.Operator(ham::Hamiltonian) = Operator(ham.basis_l, ham.data)
system(::DataOperator) = nothing
system(ham::Hamiltonian) = ham.sys
Base.:(*)(op::Operator{B1, B2}, ham::Hamiltonian{Sys, B2}) where{Sys, B1, B2} = op * Operator(ham)
Base.:(*)(ham::Hamiltonian{Sys, B2}, op::Operator{B1, B2}) where{Sys, B1, B2} = Operator(ham) * op
Base.:(+)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where{Sys, B} = op + Operator(ham)
Base.:(+)(ham::Hamiltonian{Sys, B}, op::Operator{B, B}) where{Sys, B} = Operator(ham) + op

Hamiltonian(opb::OperatorBuilder; kw...) = Hamiltonian(opb.sys, Operator(opb; kw...))

tightbinding_hamiltonian(sys::System; t1=1, t2=0, t3=0, field=NoField()) =
    build_hamiltonian(sys,
    t1 => default_bonds(sys.sample),
    t2 => default_bonds(sys.sample, Val(2)),
    t3 => default_bonds(sys.sample, Val(3)),
    field=field)
@accepts_system tightbinding_hamiltonian

const AbstractSiteOffset = Union{SiteOffset, SingleBond}
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:Tuple}; kw...)
    # Several entries
    opdata, bonds = arg
    for bond in bonds
        build_operator!(builder, sample, opdata => bond; kw...)
    end
end
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:AbstractSiteOffset};
        field=NoField())
    # Hopping operator
    opdata, bond = arg
    add_hoppings!(builder, sample.latt, opdata, bond, field)
end
function build_operator!(builder::SparseMatrixBuilder, sample::Sample, arg::Pair{<:Any, <:LatticeValue}; kw...)
    # Diagonal operator
    opdata, lv = arg
    add_diagonal!(builder, opdata, lv.values)
end
function build_operator!(builder::SparseMatrixBuilder, ::Sample, arg::DataOperator; kw...)
    # Arbitrary sparse operator
    increment!(builder, arg.data)
end

function preprocess_argument(sample::Sample, arg::DataOperator)
    if samebases(basis(arg), basis(sample))
        sparse(arg)
    elseif samebases(basis(arg), sample.internal)
        sparse(arg) ⊗ one(LatticeBasis(sample.latt))
    elseif samebases(basis(arg), LatticeBasis(sample.latt))
        if hasinternal(sample)
            sparse(arg)
        else
            internal_one(sample) ⊗ sparse(arg)
        end
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
        check_samelattice(on_lattice, sample)
        return opdata => on_lattice
    elseif on_lattice isa Number
        return preprocess_argument(sample, opdata * on_lattice)
    elseif on_lattice isa AbstractSiteOffset || on_lattice isa Tuple
        return opdata => on_lattice
    else
        error("Invalid Pair argument: unsupported on-lattice operator type")
    end
end

function build_operator(sys::System, args...; field=NoField())
    sample = sys.sample
    builder = SparseMatrixBuilder{ComplexF64}(length(sample), length(sample))
    for arg in args
        build_operator!(builder, sample, preprocess_argument(sample, arg);
            field=field)
    end
    op = Operator(basis(sample), to_matrix(builder))
    return manybodyoperator(sys, op)
end
@accepts_system build_operator
build_hamiltonian(sys::System, args...; field=NoField()) = Hamiltonian(sys, build_operator(sys, args...; field=field))
@accepts_system build_hamiltonian
