using StaticArrays

abstract type LatticeTransform end
lattransform(ltr::LatticeTransform, any) =
    throw(ArgumentError("Cannot transform object of type $(typeof(any)) with $ltr"))
lattransform(::LatticeTransform, T::Type) = T

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

struct CompositeTransform{TupleT} <: LatticeTransform
    transforms::TupleT
    CompositeTransform(transforms::LatticeTransform...) = new{typeof(transforms)}(transforms)
end
lattransform(c::CompositeTransform, a::Union{Number,SVector,SMatrix}) =
    foldl((v, t) -> lattransform(t, v), c.transforms, init=a)
Base.:(*)(lt1::LatticeTransform, lt2::LatticeTransform) = CompositeTransform(lt1, lt2)
Base.:(|>)(any, lt::LatticeTransform) = lattransform(lt, any)

lattransform(ltr::LatticeTransform, nt::NamedTuple) =
    NamedTuple{keys(nt)}(lattransform.(Ref(ltr), values(nt)))
function lattransform(ltr::LatticeTransform, lw::LatticeWithParams)
    new_lat = lw.lat |> ltr
    new_params = lw.params |> ltr
    return LatticeWithParams(new_lat, new_params)
end
lattransform(ltr::LatticeTransform, lwr::LatticeValueWrapper) =
    LatticeValueWrapper(lwr.latt |> ltr, lwr.values)
lattransform(::LatticeTransform, ::UndefinedLattice) = UndefinedLattice()
lattransform(ltr::LatticeTransform, bonds::AbstractBonds) = adapt_bonds(bonds, lattice(bonds) |> ltr)
lattransform(ltr::LatticeTransform, tr::Translation) = adapt_bonds(tr.R |> ltr, tr.lat |> ltr)
lattransform(ltr::LatticeTransform, sd::SiteDistance) = SiteDistance(sd.f, sd.lat |> ltr)
lattransform(s::Scale, sd::SiteDistance) = SiteDistance(dist -> sd.f(dist / s.scale_factor), sd.lat |> s)
lattransform(ltr::LatticeTransform, b::TwistedBoundary) = TwistedBoundary(b.translation |> ltr, b.Θ)
lattransform(ltr::LatticeTransform, bcs::BoundaryConditions) =
    BoundaryConditions(lattransform.(Ref(ltr), bcs.bcs)..., depth=bcs.depth)
lattransform(ltr::LatticeTransform, def::DefaultTranslations) =
    DefaultTranslations(def.translations |> ltr)
