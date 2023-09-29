using LinearAlgebra, Logging, StaticArrays

struct BravaisLattice{N, B<:UnitCell{Sym, N} where Sym, BS<:BoundaryConditions} <: AbstractLattice{BravaisSite{N, B}}
    bravais::B
    pointers::Vector{BravaisPointer{N}}
    boundaries::BS
    function BravaisLattice(bravais::B, pointers::Vector{BravaisPointer{N}}, boundaries::BS) where {N,B,BS}
        new{N,B,BS}(bravais, sort!(pointers), boundaries)
    end
end

BravaisLattice(bravais, pointers) = BravaisLattice(bravais, pointers, BoundaryConditions())
add_boundaries(l::BravaisLattice, bs) = BravaisLattice(l.bravais, l.pointers, to_boundaries(bs))
sites(l::BravaisLattice) = Sites(add_boundaries(l, BoundaryConditions()))

Base.:(==)(l1::BravaisLattice, l2::BravaisLattice) = (l1.pointers == l2.pointers) && (l1.bravais == l2.bravais)

Base.emptymutable(l::BravaisLattice{B, N}, ::Type{BravaisSite{N}}=eltype(l)) where {B, N} =
    BravaisLattice(l.bravais, [])
Base.copymutable(l::BravaisLattice) = BravaisLattice(l.bravais, copy(l.pointers))
Base.length(l::BravaisLattice) = length(l.pointers)
lattice_type(::BravaisLattice{<:UnitCell{Sym}}) where {Sym} = Sym
basis_length(l::BravaisLattice) = length(l.bravais)
check_samebravais(l::BravaisLattice, site::BravaisSite) =
    @assert l.bravais == site.bravais

default_bonds(::BravaisLattice, ::Val) = ()
default_bonds(l::BravaisLattice) = default_bonds(l, Val(1))
default_bonds(l::BravaisLattice, i::Int) = default_bonds(l, Val(i))

get_site(l::BravaisLattice, lp::BravaisPointer) = BravaisSite(lp, l.bravais)
get_site(::BravaisLattice, site::BravaisSite) = site
get_site(::BravaisLattice, ::NoSite) = NoSite()

shift_site(l::BravaisLattice, site::AbstractSite) = shift_site(l.boundaries, l, site)

Base.in(lp::BravaisPointer, l::BravaisLattice) = insorted(lp, l.pointers)
function Base.in(site::BravaisSite{N, B}, l::BravaisLattice{N, B}) where {N, B}
    @boundscheck check_samebravais(l, site)
    return in(site.lp, l)
end

Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, i::Int)
    @boundscheck checkbounds(l, i)
    return get_site(l, l.pointers[i])
end
Base.@propagate_inbounds function Base.getindex(l::BravaisLattice, is::AbstractVector{Int})
    @boundscheck checkbounds(l, is)
    return BravaisLattice(l.bravais, l.pointers[sort(is)])
end
Base.@propagate_inbounds function Base.delete!(l::BravaisLattice, lp::BravaisPointer)
    i = site_index(l, lp)
    i !== nothing && deleteat!(l.pointers, i)
    return l
end
Base.@propagate_inbounds function Base.deleteat!(l::BravaisLattice, inds)
    @boundscheck checkbounds(l, inds)
    deleteat!(l.pointers, inds)
    return l
end

function Base.push!(l::BravaisLattice, lp::BravaisPointer)
    @assert 1 ≤ lp.basis_index ≤ basis_length(lp) "invalid basis index $(lp.basis_index)"
    i = searchsortedfirst(l.pointers, lp)
    i < length(l) && l.pointers[i] != lp && insert(l.pointers, i, lp)
end
Base.push!(l::BravaisLattice{N, B}, site::BravaisSite{N, B}) where {N, B} =
    push!(l, site.lp)

Base.@propagate_inbounds function site_index(l::BravaisLattice, lp::BravaisPointer)
    i = searchsortedfirst(l.pointers, lp)
    i > length(l) && return nothing
    return l.pointers[i] == lp ? i : nothing
end
Base.@propagate_inbounds function site_index(l::BravaisLattice, site::BravaisSite)
    @boundscheck check_samebravais(l, site)
    site_index(l, site.lp)
end

function Base.iterate(l::BravaisLattice, state = (1, length(l)))
    i, len = state
    return i > len ? nothing : (l[i], (i+1, len))
end

function Base.show(io::IO, ::MIME"text/plain", l::BravaisLattice{N, <:UnitCell{Sym}}) where {N,Sym}
    print(io, length(l), "-site ", N, "-dimensional ", Sym, " lattice")
    if basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

function BravaisLattice(uc::UnitCell{Sym,N,NB}, sz::NTuple{N,Int};
        boundaries=BoundaryConditions(),
        offset = @SVector zeros(N)) where {Sym,N,NB}
    ptrs = BravaisPointer{N}[]
    for unit_cell in CartesianIndex{N}():CartesianIndex(reverse(sz))
        svec = reverse(SVector{N}(Tuple(unit_cell)))
        for i in 1:NB
            push!(ptrs, BravaisPointer(svec, i))
        end
    end
    if offset === :center
        basis_com = vec(sum(uc.basis, dims=2) / NB)
        center_offset = uc.translation_vectors * SVector{N}(@. (sz + 1) / 2) + basis_com
        BravaisLattice(UnitCell(uc, -center_offset), ptrs, to_boundaries(boundaries))
    elseif offset isa AbstractVector{<:Number}
        BravaisLattice(UnitCell(uc, offset), ptrs, to_boundaries(boundaries))
    else
        throw(ArgumentError("invalid `offset` keyword argument"))
    end
end

const InfDimLattice{Sym,N,NB} = BravaisLattice{N, <:UnitCell{Sym,N,NB}}
const AnyDimLattice{Sym,NB} = BravaisLattice{N, <:UnitCell{Sym,N,NB}} where N
function InfDimLattice{Sym,N,NB}(sz::Vararg{Int,N}; kw...) where {Sym,N,NB}
    return BravaisLattice(UnitCell{Sym,N,NB}(), sz; kw...)
end
function AnyDimLattice{Sym,NB}(sz::Vararg{Int,N}; kw...) where {Sym,N,NB}
    return BravaisLattice(UnitCell{Sym,N,NB}(), sz; kw...)
end
function (::Type{T})(f::Function, args...; kw...) where T<:BravaisLattice
    sublattice(f, T(args...; kw...))
end

macrocell_size(::BravaisLattice) = error("This function is discontinued")
