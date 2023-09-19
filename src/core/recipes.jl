using RecipesBase

@recipe function f(site::AbstractSite)
    seriestype := :scatter
    [Tuple(site.coords),]
end

@recipe function f(l::AbstractLattice, style::Symbol=:plain)
    l, Val(style)
end

@recipe function f(::AbstractLattice, ::Val{StyleT}) where StyleT
    error("Unsupported lattice plot style $StyleT")
end

@recipe function f(l::AbstractLattice, ::Val{:plain})
    l, nothing
end

@recipe function f(l::AbstractLattice, ::Val{:high_contrast})
    markersize := 4
    markercolor := :black
    markerstrokealpha := 1
    markerstrokestyle := :solid
    markerstrokewidth := 2
    markerstrokecolor := :white
    l, nothing
end

@recipe function f(l::AbstractLattice, ::Val{:markers})
    label --> ""
    @series begin   # The sites
        seriestype := :scatter
        annotations = [(" " * string(i), :left, :top, :grey, 6) for i in eachindex(l)]
        series_annotations := annotations
        l, nothing
    end
end

@recipe function f(l::AbstractLattice, ::Val{:pretty})
    @series l, Val(:markers)
    @series l, Val(:bonds)
end

@recipe function f(l::AbstractLattice, v)
    aspect_ratio := :equal
    marker_z := v
    pts = collect_coords(l)
    if v !== nothing && RecipesBase.is_key_supported(:hover)
        crd_markers = [join(round.(crds[:, i], digits=3), ", ") for i in 1:size(pts, 2)]
        hover := string.(round.(v, digits=3), " @ (", crd_markers, ")")
    end
    if dims(l) == 1
        seriestype := :scatter
        @series vec(pts), zeros(vec(pts))
    elseif dims(l) == 2
        seriestype --> :scatter
        X, Y = eachrow(pts)
        if plotattributes[:seriestype] == :scatter
            @series X, Y
        elseif plotattributes[:seriestype] == :surface
            @series X, Y, v
        else
            throw(ArgumentError("series type $(plotattributes[:seriestype]) not supported for 2D lattices"))
        end
    elseif dims == 3
        seriestype := :scatter3d
        @series Tuple(eachrow(pts))
    end
end

@recipe function f(lv::LatticeValue{<:Number})
    lv.lattice, lv.values
end

@recipe function f(ag::AdjacencyMatrix)
    aspect_ratio := :equal
    l = lattice(ag)
    pts = NTuple{dims(l), Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for i in eachindex(l), j in eachindex(l)
        site1 = l[i]
        site2 = l[j]
        !ag[site1, site2] && continue
        A = site1.coords
        B = site2.coords
        push!(pts, Tuple(A), Tuple(B), br_pt)
    end
    label := nothing
    pts
end
