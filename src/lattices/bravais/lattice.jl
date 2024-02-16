using LinearAlgebra, Logging, StaticArrays

struct BravaisLattice{N, B, BS<:BoundaryConditions} <: AbstractLattice{BravaisSite{N, B}}
    bravais::B
    pointers::Vector{BravaisPointer{N}}
    boundaries::BS
    b_depth::Int
    function BravaisLattice(bravais::B, pointers::Vector{BravaisPointer{N}}, boundaries::BS; b_depth=1) where {N,B,BS}
        dims(bravais) != N && throw(ArgumentError("Dimension mismatch"))
        new{N,B,BS}(bravais, sort!(pointers), boundaries, b_depth)
    end
end
BravaisLattice(bravais, pointers) = BravaisLattice(bravais, pointers, BoundaryConditions())
add_boundaries(l::BravaisLattice, bs) = BravaisLattice(l.bravais, l.pointers, to_boundaries(bs))
add_boundaries(l::BravaisLattice, ::Nothing) = l
add_boundaries(l::AbstractLattice, ::Nothing) = l
sites(l::BravaisLattice) = Sites(add_boundaries(l, BoundaryConditions()))
sublatvector(l::BravaisLattice, i::Int) = sublatvector(l.bravais, i)
trvectors(l::BravaisLattice) = trvectors(l.bravais)
b_depth(l::BravaisLattice) = l.b_depth
cartesian_indices(l::BravaisLattice{N, B, <:BoundaryConditions{<:NTuple{M}}} where {N, B}) where M =
    cartesian_indices(l.b_depth, Val(M))

Base.:(==)(l1::BravaisLattice, l2::BravaisLattice) = (l1.pointers == l2.pointers) && (l1.bravais == l2.bravais)

Base.emptymutable(l::BravaisLattice{N, B}, ::Type{BravaisSite{N, B}}) where {N, B} =
    BravaisLattice(l.bravais, BravaisPointer{N}[])
Base.copymutable(l::BravaisLattice) = BravaisLattice(l.bravais, copy(l.pointers))
Base.length(l::BravaisLattice) = length(l.pointers)
lattice_type(::BravaisLattice{<:UnitCell{Sym}}) where {Sym} = Sym
basis_length(l::BravaisLattice) = length(l.bravais)
check_samebravais(l::BravaisLattice, site::BravaisSite) =
    @assert l.bravais == site.bravais

get_site(l::BravaisLattice, lp::BravaisPointer) = BravaisSite(lp, l.bravais)
get_site(::BravaisLattice, site::BravaisSite) = site
get_site(::BravaisLattice, ::NoSite) = NoSite()

function resolve_site(l::BravaisLattice, site::AbstractSite)
    factor, site = shift_site(l.boundaries, l, site)
    i = site_index(l, site)
    i === nothing && return nothing
    return ResolvedSite(site, i, factor)
end

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
    @assert 1 ≤ lp.basis_index ≤ basis_length(l) "invalid basis index $(lp.basis_index)"
    i = searchsortedfirst(l.pointers, lp)
    i ≤ length(l) && l.pointers[i] == lp && return l
    insert!(l.pointers, i, lp)
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

function Base.summary(io::IO, l::BravaisLattice{N, <:UnitCell{Sym}}) where {N,Sym}
    print(io, length(l), "-site ", N, "-dim ", Sym, " lattice")
    if basis_length(l) > 1
        print(io, " (", basis_length(l), "-site basis)")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", l::BravaisLattice)
    summary(io, l)
    (isempty(l.boundaries.bcs) || get(io, :compact, false)) && return
    print(io, " with boundary conditions:")
    for bc in l.boundaries.bcs
        println()
        show(io, mime, bc)
    end
end

function span_unitcells(uc::UnitCell{Sym,N,NB}, sz::NTuple{N,Int};
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
span_unitcells(uc::UnitCell, sz::Tuple{Vararg{Int}}; kw...) =
    throw(ArgumentError("Dimension mismatch: $(dims(uc))-dim unit cell, $(length(sz)) lattice dimensions"))
