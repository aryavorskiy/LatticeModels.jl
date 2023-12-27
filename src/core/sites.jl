struct Sites{SiteT,LT<:AbstractLattice{SiteT}}<:AbstractLattice{SiteT}
    latt::LT
    Sites(l::LT) where {SiteT,LT<:AbstractLattice{SiteT}} = new{SiteT, LT}(l)
end
sites(s::Sites) = s
sites(l::AbstractLattice) = Sites(l)
sites(any) = sites(lattice(any))
Base.:(==)(s1::Sites, s2::Sites) = s1.latt == s2.latt
Base.:(==)(s::Sites, l::AbstractLattice) = s.latt == l
Base.:(==)(l::AbstractLattice, s::Sites) = s.latt == l
Base.length(s::Sites) = length(s.latt)
Base.push!(s::Sites, arg) = (Base.push!(s.latt, arg); s)
Base.iterate(s::Sites, state...) = iterate(s.latt, state...)
Base.@propagate_inbounds Base.getindex(s::Sites, any...; kw...) = getindex(s.latt, any...; kw...)
Base.@propagate_inbounds site_index(s::Sites, any) = site_index(s.latt, any)
Base.deleteat!(s::Sites, i) = (Base.deleteat!(s.latt, i); s)

function Base.summary(io::IO, sites::Sites)
    print(io, "Sites(")
    summary(io, sites.latt)
    print(io, ")")
end
