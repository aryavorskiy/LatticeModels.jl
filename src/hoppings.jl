import Base: show, |

struct Hopping{MT<:AbstractMatrix}
    site_indices::Tuple{Int, Int}
    tr_vector::Vector{Int}
    pbc::Vector{Bool}
    hop_operator::MT
    function Hopping(site_indices, tr_vector, pbc, hop_operator)
        @assert length(tr_vector) == length(pbc) "inconsistent dimensionality"
        new{typeof(hop_operator)}(site_indices, tr_vector, pbc, hop_operator)
    end
end

_wrap_operator(o::Number) = [o;;]
_wrap_operator(o::AbstractMatrix) = o
_tr_and_pbc(tr_vc::Vector, pbc::Vector) = (tr_vc, pbc)
_tr_and_pbc(tr_vc::Vector, pbc::Bool) = (tr_vc, fill(pbc, length(tr_vc)))
_tr_and_pbc(pbc::Vector) = (zeros(Int, length(pbc)), pbc)
_tr_and_pbc(pbc::Bool) = ([0], [pbc])
function _tr_and_pbc(axis::Int, pbc::Vector)
    tr_vc = zeros(Int, length(pbc))
    tr_vc[axis] = 1
    return (tr_vc, pbc)
end
function _tr_and_pbc(axis::Int, pbc::Bool)
    tr_vc = zeros(Int, axis)
    tr_vc[axis] = 1
    return (tr_vc, fill(pbc, axis))
end
function Hopping(hop_operator=1;site_indices::NTuple{2} = (1, 1), pbc=false, kw...)
    if :axis in keys(kw)
        tr_vc, pbc = _tr_and_pbc(kw[:axis], pbc)
    elseif :tr_vector in keys(kw)
        tr_vc, pbc = _tr_and_pbc(kw[:tr_vector], pbc)
    elseif site_indices[1] != site_indices[2]
        tr_vc, pbc = _tr_and_pbc(pbc)
    else throw(ArgumentError("hopping connects site to itself"))
    end
    Hopping(site_indices, tr_vc, pbc, _wrap_operator(hop_operator))
end

function show(io::IO, m::MIME"text/plain", hop::Hopping)
    println(io, "Hopping")
    println(io, "Connects site #$(hop.site_indices[1]) with site #$(hop.site_indices[1]) translated by $(hop.tr_vector)")
    if !iszero(hop.tr_vector)
        print(io, "Boundary conditions: ")
        if all(hop.pbc .| (hop.tr_vector .== 0))
            println(io, "periodic")
        else
            p_axes = [i for i in eachindex(hop.tr_vector) if hop.pbc[i] || hop.tr_vector[i] == 0]
            if length(p_axes) == 0
                println(io, "open")
            elseif length(p_axes) == 1
                println(io, "periodic at axis $(only(p_axes))")
            else
                println(io, "periodic at axes $(join(p_axes, "×"))")
            end
        end
    end
    print(io, "Hopping operator matrix: ")
    show(io, m, hop.hop_operator)
end
dims(h::Hopping) = length(h.tr_vector)
dims_internal(h::Hopping) = size(h.hop_operator)[1]

function _promote_dims!(h::Hopping, ndims::Int)
    if ndims ≥ dims(h)
        append!(h.pbc, fill(false, ndims - dims(h)))
        append!(h.tr_vector, fill(0, ndims - dims(h)))
    else
        for _ in 1:dims(h) - ndims
            if h.tr_vector[end] == 0
                pop!(h.tr_vector)
                pop!(h.pbc)
            else
                throw(ArgumentError("cannot shrink hopping dims, non-zero translation found"))
            end
        end
    end
end

function _match(h::Hopping, l::AbstractLattice, site1::LatticeIndex, site2::LatticeIndex)
    ret = (site1.basis_index, site2.basis_index) == h.site_indices
    for i in 1:dims(h)
        vi = site2.unit_cell[i] - site1.unit_cell[i] - h.tr_vector[i]
        if h.pbc[i] && (vi % _sz(l)[i] != 0)
            ret = false
        elseif !h.pbc[i] && (vi != 0)
            ret = false
        end
    end
    return ret
end

function _hopping_operator!(lop::LatticeOperator, lf::Function, l::AbstractLattice, hop::Hopping, field::AbstractField)
    _promote_dims!(hop, dims(l))
    d = dims(l)
    p1 = zeros(d)
    p2 = zeros(d)
    buf = zeros(d)
    bvs = bravais(l)
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            if @inbounds(_match(hop, l, site1, site2)) && lf(l, site1)
                coords!(p1, l, site1, bvs, buf)
                copy!(p2, p1)
                p2 .+= hop.tr_vector
                pmod = exp(2π * im * _trip_integral!(field, p1, p2, buf))
                lop[i, j] += hop.hop_operator * pmod
                lop[j, i] += hop.hop_operator' * pmod'
            end
            j += 1
        end
        i += 1
    end
    lop
end

function hopping_operator(f::Function, l::AbstractLattice, hop::Hopping; field::AbstractField=NoField())
    lop = _zero_on_basis(l, hop.hop_operator)
    _hopping_operator!(lop, _propagate_lattice_args(f, l), l, hop, field)
    return lop
end

function hopping_operator(l::AbstractLattice, hop::Hopping; field::AbstractField=NoField())
    lop = _zero_on_basis(l, hop.hop_operator)
    _hopping_operator!(lop, _always_true_on_lattice, l, hop, field)
    return lop
end

macro hopping_operator(for_loop::Expr)
    if for_loop.head !== :for
        throw(ArgumentError("expression must be a for loop, not $(for_loop.head)"))
    elseif length(for_loop.args) != 2
        throw(ArgumentError("malformed for loop")) # This should never happen, but still...
    end
    itr::Expr, body::Expr = for_loop.args
    if !Meta.isexpr(itr, :(=), 2)
        throw(ArgumentError("invalid for loop iteration specification; must be a simple assignment"))
    end
    itr_vars, lattice_var = itr.args
    if !Meta.isexpr(itr_vars, :tuple, 2)
        throw(ArgumentError("invalid for loop iterator variable; must be a 2-tuple"))
    end
    site1_var, site2_var = itr_vars.args
    while body.args[end] isa LineNumberNode
        pop!(body.args)
    end
    dump(body)
    quote
        i = 0
        l = $(esc(lattice_var))
        local matrix = nothing
        local N = 0
        for $(esc(site1_var)) in l
            i += 1
            j = 0
            for $(esc(site2_var)) in l
                j += 1
                if i ≥ j
                    continue
                end
                block_res = $(esc(body))
                if block_res !== nothing
                    if matrix === nothing
                        if size(block_res) == ()
                            N = 1
                        else
                            N = size(block_res)[1]
                        end
                        matrix = zeros(ComplexF64, N * length(l), N * length(l))
                    end
                    matrix[N * (i - 1) + 1: N * i, N * (j - 1) + 1: N * j] .= block_res
                    matrix[N * (j - 1) + 1: N * j, N * (i - 1) + 1: N * i] .= block_res'
                end
            end
        end
        LatticeVecOrMat(Basis(l, N), matrix)
    end
end

struct BondSet
    lattice::AbstractLattice
    bmat::BitMatrix
    function BondSet(l::AbstractLattice, bmat::BitMatrix)
        @assert all(size(bmat) .== length(l)) "inconsistent connectivity matrix size"
        new(l, bmat)
    end
    BondSet(l::AbstractLattice, bmat::AbstractMatrix{Bool}) = new(l, convert(BitMatrix, bmat))
    function BondSet(l::AbstractLattice)
        new(l, zeros(Bool, length(l), length(l)))
    end
end

function BondSet(op::LatticeOperator)
    N = dims_internal(op)
    matrix = Bool[!iszero(op[i, j])
                    for i in 1:length(op.basis.lattice), j in 1:length(op.basis.lattice)]
    return BondSet(op.basis.lattice, matrix)
end

function BondSet(l::AbstractLattice, hops::Hopping...)
    matrix = zeros(Bool, length(l), length(l))
    for h in hops
        _promote_dims!(h, dims(l))
    end
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            matrix[i, j] = any(@inbounds(_match(h, l, site1, site2)) for h in hops)
            j += 1
        end
        i += 1
    end
    return BondSet(l, matrix)
end

function |(bss::BondSet...)
    @assert allequal(getproperty.(bss, :lattice))
    return BondSet(bss[1].lattice, .|(getproperty.(bss, :bmat)...))
end

function show(io::IO, m::MIME"text/plain", bs::BondSet)
    println(io, "BondSet with $(count(==(true), bs.bmat)) bonds")
    print(io, "on ")
    show(io, m, bs.lattice)
end

@recipe function f(bs::BondSet)
    aspect_ratio := :equal
    l = bs.lattice
    pts = Tuple{Float64, Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    i = 1
    for site1 in l
        j = 1
        A = coords(l, site1)
        for site2 in l
            if i != j && bs.bmat[i, j]
                B = coords(l, site2)
                T = radius_vector(l, site2, site1)
                push!(pts, Tuple(A))
                push!(pts, Tuple(A + T/2))
                push!(pts, br_pt)
                push!(pts, Tuple(B))
                push!(pts, Tuple(B - T/2))
                push!(pts, br_pt)
            end
            j += 1
        end
        i += 1
    end
    @series begin
        label := nothing
        pts
    end
end
