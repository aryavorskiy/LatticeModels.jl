using RecipesBase, LinearAlgebra, Logging
import Base: length, size, iterate, show

struct Bravais
    translation_vectors::Matrix{Float64}
    basis::Matrix{Float64}
    function Bravais(translation_vectors::Matrix{<:Real}, basis::Matrix{<:Real})
        @assert size(basis)[1] == size(translation_vectors)[1]
        new(translation_vectors .|> Float64, basis .|> Float64)
    end
    function Bravais(translation_vectors::Matrix{<:Real})
        new(translation_vectors, zeros(Float64, (size(translation_vectors)[1], 1)))
    end
end

dims(b::Bravais) = size(b.basis)[1]
length(b::Bravais) = size(b.basis)[2]

abstract type AbstractLattice end
abstract type FiniteBravaisLattice<:AbstractLattice end

size(l::FiniteBravaisLattice) = l.lattice_size
dims(l::FiniteBravaisLattice) = length(l.lattice_size)
bravais(l::LT) where LT<:FiniteBravaisLattice = _bravais(l, dims(l))

length(l::FiniteBravaisLattice) = prod(size(l)) * length(bravais(l))
mutable struct LatticeIndex
    unit_cell::Vector{Int}
    basis_index::Int
end

start(l::FiniteBravaisLattice) = LatticeIndex(fill(1, dims(l)), 1)

function proceed!(l::FiniteBravaisLattice, bvs_len, site::LatticeIndex)
    @inbounds if site.basis_index == bvs_len
        site.basis_index = 1
        i = length(site.unit_cell)
        while site.unit_cell[i] == size(l)[i]
            site.unit_cell[i] = 1
            i -= 1
            i == 0 && return false
        end
        site.unit_cell[i] += 1
    else
        site.basis_index += 1
    end
    return true
end

function iterate(l::FiniteBravaisLattice)
    site = start(l)
    return (site, (site, length(bravais(l)), 1))
end

function iterate(l::FiniteBravaisLattice, state::Tuple{LatticeIndex, Int, Int})
    site, bvs_len, index = state
    !proceed!(l, bvs_len, site) && return nothing
    index += 1
    return (site, (site, bvs_len, index))
end

function coords!(crd::Vector, l::AbstractLattice, site::LatticeIndex, bvs::Bravais, buf::Vector)
    copy!(buf, site.unit_cell)
    buf .-= _sz(l) ./ 2
    buf .-= 0.5
    mul!(crd, bvs.translation_vectors, buf)
    crd += bvs.basis[:, site.basis_index]
end

coords(l::AbstractLattice, site::LatticeIndex) =
    bravais(l).basis[:, site.basis_index] +  bravais(l).translation_vectors * (site.unit_cell - _sz(l) / 2 .- 0.5)

_bravais(::FiniteBravaisLattice, ::Int) = throw(ArgumentError("No Bravais lattice defined for this type"))
_assert_dims(l::LT, ::Int) where LT<:FiniteBravaisLattice = l
_prettyprint_name(::LT) where {LT<:FiniteBravaisLattice} = string(LT)

function _propagate_lattice_args(f, l::AbstractLattice)
    if hasmethod(f, NTuple{dims(l), Number})
        (_l::AbstractLattice, site::LatticeIndex) -> f(coords(_l, site)...)
    elseif hasmethod(f, Tuple{LatticeIndex})
        (::AbstractLattice, site::LatticeIndex) -> f(site)
    elseif hasmethod(f, Tuple{LatticeIndex, Vararg{Number, dims(l)}})
        (_l::AbstractLattice, site::LatticeIndex) -> f(site, coords(_l, site)...)
    else
        throw(ArgumentError("failed to propagate args: unsupported lambda type"))
    end
end

_always_true_on_lattice(::AbstractLattice, ::LatticeIndex) = true

macro lattice_def(struct_block::Expr)
    if struct_block.head != :struct
        error("Struct block expected")
    end
    struct_name = struct_block.args[2]
    if !(struct_name isa Symbol)
        error("Invalid struct name")
    end
    body = struct_block.args[3]
    struct_definition = quote
        struct $struct_name <: FiniteBravaisLattice
            lattice_size::Vector{Int}
            function $struct_name(lattice_size)
                _assert_dims(new(lattice_size |> collect), length(lattice_size))
            end
            $struct_name(lattice_sz_by_axes::Int...) = $struct_name(lattice_sz_by_axes)
        end
        $struct_name(f::Function, args...) = sublattice(f, $struct_name(args...))
        export $struct_name
    end
    bravais_flag = false
    for expr in body.args
        !(expr isa Expr) && continue
        if expr.head == :(:=)
            param, value = expr.args
            if param == :name
                if value isa String
                    push!(struct_definition.args, :(_prettyprint_name(::$struct_name) = $value))
                else
                    @warn "lattice name should be a string; ignored"
                end
            end
        elseif bravais_flag
            @warn "multiple Bravais lattice definitions detected; ignored"
            continue
        elseif expr.head == :-> && expr.args[1] isa Symbol
            fn_arg, fn_body = expr.args
            push!(struct_definition.args, :(_bravais(::$struct_name, $fn_arg::Int) = $fn_body))
            bravais_flag = true
        elseif expr.head == :call && expr.args[1] == :Bravais
            call_args = expr.args[2]
            push!(struct_definition.args, :(bravais(::$struct_name) = $expr))
            if call_args.head == :vcat
                message = "lattice type $(string(struct_name)) only supports dimension count $(length(call_args.args))"
                push!(struct_definition.args, :(
                function _assert_dims(l::$struct_name, dims::Int)
                    @assert dims == $(length(call_args.args)) $message
                    return l
                end))
            end
            bravais_flag = true
        elseif expr.head == :vcat
            push!(struct_definition.args, :(bravais(::$struct_name) = Bravais($expr)))
            message = "lattice type $(string(struct_name)) only supports dimension count $(length(expr.args))"
            push!(struct_definition.args, :(
            function _assert_dims(l::$struct_name, dims::Int)
                @assert dims == $(length(expr.args)) $message
                return l
            end))
            bravais_flag = true
        end
    end
    return :($(esc(struct_definition)))
end

@lattice_def struct SquareLattice
    name := "square lattice"
    dims -> Bravais(Matrix(I, (dims, dims)))
end

@lattice_def struct TriangularLattice
    name := "triangular lattice"
    [1 0.5;0 √3/2]
end

@lattice_def struct HoneycombLattice
    name := "honeycomb lattice"
    Bravais([1 0.5;0 √3/2], [0 0.5;0 √3/6])
end

struct SubLattice{LT}<:AbstractLattice where LT<:FiniteBravaisLattice
    lattice::LT
    mask::Vector{Bool}
    function SubLattice(lattice::T, mask::Vector{Bool}) where T
        @assert length(mask) == length(lattice) "Lattice site count does not match mask length"
        new{T}(lattice, mask)
    end
end

function sublattice(f::Function, l::AbstractLattice)
    lf = _propagate_lattice_args(f, l)
    SubLattice(l, [lf(l, site) for site in l])
end

dims(sl::SubLattice) = dims(sl.lattice)
bravais(sl::SubLattice) = bravais(sl.lattice)
length(sl::SubLattice) = count(sl.mask)
coords(sl::SubLattice, state::LatticeIndex) = coords(sl.lattice, state)

function iterate(sl::SubLattice)
    site = start(sl.lattice)
    bvs_len = length(bravais(sl))
    index = 1
    while !sl.mask[index]
        proceed!(sl.lattice, bvs_len, site)
        index += 1
    end
    return (site, (site, bvs_len, index))
end

function iterate(sl::SubLattice, state::Tuple{LatticeIndex, Int, Int})
    site, bvs_len, index = state
    !proceed!(sl.lattice, bvs_len, site) && return nothing
    index += 1
    while !sl.mask[index]
        !proceed!(sl.lattice, bvs_len, site) && return nothing
        index += 1
    end
    return (site, (site, bvs_len, index))
end

_sz(l::FiniteBravaisLattice) = size(l)
_sz(l::SubLattice) = size(l.lattice)

function radius_vector(l::AbstractLattice, site1::LatticeIndex, site2::LatticeIndex)
    ret_vec = coords(l, site1) - coords(l, site2)
    tr_diff = (site1.unit_cell - site2.unit_cell) ./ _sz(l)
    bravais_tr_vecs = bravais(l).translation_vectors
    for i in eachindex(tr_diff)
        if tr_diff[i] > 0.5
            ret_vec -= bravais_tr_vecs[:, i] * _sz(l)[i]
        elseif tr_diff[i] < -0.5
            ret_vec += bravais_tr_vecs[:, i] * _sz(l)[i]
        end
    end
    return ret_vec
end

@recipe function f(sl::SubLattice)
    show_excluded --> true
    if plotattributes[:show_excluded]
        opacity := sl.mask * 0.9 .+ 0.1
        sl.lattice
    else
        sl, nothing
    end
end

@recipe function f(l::AbstractLattice, v=nothing)
    label --> nothing
    aspect_ratio := :equal
    marker_z := v
    markerstrokewidth := 0
    d = dims(l)
    pts = zeros(Float64, d, length(l))
    bvs = bravais(l)
    crd = zeros(d)
    buf = zeros(d)
    i = 1

    for site in l
        coords!(crd, l, site, bvs, buf)
        pts[:, i] = crd
        i += 1
    end
    if dims(l) == 1
        pts = vcat(pts, zeros(Float64, 1, length(l)))
    end
    if dims(l) == 3
        X, Y, Z = eachrow(pts)
        Xr, Yr, Zr = eachrow(round.(pts, digits=3))
        seriestype := :scatter3d
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ", ", Zr, ")")
        end
        X, Y, Z
    else
        X, Y = eachrow(pts[1:2, :])
        Xr, Yr = eachrow(round.(pts[1:2, :], digits=3))
        seriestype --> :scatter
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ")")
        end
        if plotattributes[:seriestype] == :scatter
            X, Y
        elseif plotattributes[:seriestype] == :surface
            X, Y, v
        else throw(ArgumentError("unsupported series type $(plotattributes[:seriestype])"))
        end
    end
end

function show(io::IO, ::MIME"text/plain", l::FiniteBravaisLattice)
    print(io, join(size(l), "×") * " " * _prettyprint_name(l))
end

function show(io::IO, m::MIME"text/plain", sl::SubLattice)
    print(io, "$(length(sl))-element SubLattice of ")
    show(io, m, sl.lattice)
end
