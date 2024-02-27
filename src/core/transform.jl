using StaticArrays

abstract type AbstractLatticeTransform end
abstract type LatticeTransform <: AbstractLatticeTransform end
lattransform(ltr::AbstractLatticeTransform, any) =
    throw(ArgumentError("Cannot transform object of type $(typeof(any)) with $ltr"))
lattransform(::AbstractLatticeTransform, T::Type) = T
struct Shift{N} <: LatticeTransform
    shift_vector::SVector{N,Float64}
    Shift(vec::AbstractVector) = new{length(vec)}(SVector{length(vec)}(vec))
end

lattransform(::Shift, n::Number) = n
lattransform(s::Shift{N}, vec::SVector{N}) where N = vec + s.shift_vector
lattransform(s::Shift{N}, mat::SMatrix{N}) where N = mat .+ s.shift_vector

struct Rotate{N} <: LatticeTransform
    rotation_matrix::SMatrix{N,N,Float64}
end
Rotate(;θ::Real) = Rotate{2}(SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ)))

lattransform(::Rotate, n::Number) = n
lattransform(r::Rotate{N}, vec::SVector{N}) where N = r.rotation_matrix * vec
lattransform(r::Rotate{N}, mat::SMatrix{N}) where N = r.rotation_matrix * mat

struct Scale <: LatticeTransform
    scale_factor::Float64
end
lattransform(s::Scale, n::Number) = s.scale_factor * n
lattransform(s::Scale, vec::SVector) = s.scale_factor * vec
lattransform(s::Scale, mat::SMatrix) = s.scale_factor * mat

struct Project{ParamT<:SiteProperty} <: LatticeTransform
    projection::ParamT
end
Project(sym::Symbol) = Project(SitePropertyAlias{sym}())
Project(ps...) = *((Project(p) for p in ps)...)
lattransform(::Project, n::Number) = n
function lattransform(p::Project{Coord}, vec::SVector)
    p.projection.axis in eachindex(vec) ||
        throw(ArgumentError("Invalid axis $(p.projection.axis) for vector of length $(length(vec))"))
    vec[p.projection.axis .!= eachindex(vec)]
end
function lattransform(p::Project{Coord}, mat::SMatrix{M, N}) where {M, N}
    p.projection.axis in 1:M ||
        throw(ArgumentError("Invalid axis $(p.projection.axis) for matrix of size $(join("×", size(mat)))"))
    SMatrix{M-1,N}(mat[p.projection.axis .!= 1:size(mat, 1), :])
end

struct CompositeTransform{TupleT} <: AbstractLatticeTransform
    transforms::TupleT
    CompositeTransform(transforms::LatticeTransform...) = new{typeof(transforms)}(transforms)
end
lattransform(c::CompositeTransform, a) =
    foldl((v, t) -> lattransform(t, v), c.transforms, init=a)
Base.:(*)(lt1::LatticeTransform, lt2::LatticeTransform...) = CompositeTransform(lt1, lt2...)
Base.:(*)(lt1::CompositeTransform, lt2::CompositeTransform) = CompositeTransform(lt1.transforms..., lt2.transforms...)
Base.:(|>)(any, lt::AbstractLatticeTransform) = lattransform(lt, any)

lattransform(ltr::LatticeTransform, nt::NamedTuple) =
    NamedTuple{keys(nt)}(lattransform.(Ref(ltr), values(nt)))
function lattransform(ltr::LatticeTransform, lw::LatticeWithParams)
    new_lat = lw.lat |> ltr
    new_params = lw.params |> ltr
    return LatticeWithParams(new_lat, new_params)
end
function lattransform(proj::Project, lw::LatticeWithParams)
    new_lat = lw.lat |> proj
    new_params = lw.params |> proj
    return LatticeWithParams(new_lat, Base.structdiff(new_params, NamedTuple{(:latticetype,)}))
end
lattransform(ltr::LatticeTransform, lwr::LatticeValueWrapper) =
    LatticeValueWrapper(lwr.latt |> ltr, lwr.values)
function lattransform(proj::Project, lwr::LatticeValueWrapper)
    new_latt = lwr.latt |> proj
    new_vals = Array{eltype(lwr.values)}(undef, length(new_latt))
    for (i, site) in enumerate(lwr.latt)
        j = site_index(new_latt, site |> proj)
        new_vals[j] = lwr.values[i]
    end
    return LatticeValueWrapper(new_latt, new_vals)
end
lattransform(::LatticeTransform, ::UndefinedLattice) = UndefinedLattice()
lattransform(ltr::LatticeTransform, bonds::AbstractBonds) = adapt_bonds(bonds, lattice(bonds) |> ltr)
lattransform(ltr::LatticeTransform, tr::Translation) = adapt_bonds(tr.R |> ltr, tr.lat |> ltr)
lattransform(ltr::LatticeTransform, sd::SiteDistance) = SiteDistance(sd.f, sd.lat |> ltr)
lattransform(s::Scale, sd::SiteDistance) = SiteDistance(dist -> sd.f(dist / s.scale_factor), sd.lat |> s)
lattransform(ltr::LatticeTransform, b::TwistedBoundary) = TwistedBoundary(b.translation |> ltr, b.Θ)
lattransform(ltr::LatticeTransform, bcs::BoundaryConditions) =
    BoundaryConditions(lattransform.(Ref(ltr), bcs.bcs)..., depth=bcs.depth)
function lattransform(ltr::LatticeTransform, def::DefaultTranslations)
    processed_pairs = Tuple(pairs(def.translations |> ltr))
    filtered_pairs = filter(p -> p[2] !== nothing, processed_pairs)
    DefaultTranslations(NamedTuple(filtered_pairs))
end
