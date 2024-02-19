using RecipesBase

@recipe function f(site::AbstractSite)
    seriestype := :scatter
    [Tuple(site.coords),]
end

@recipe function f(l::AbstractLattice, style::Symbol)
    l, Val(style)
end

@recipe function f(::AbstractLattice, ::Val{StyleT}) where StyleT
    error("Unsupported lattice plot style '$StyleT'")
end

@recipe function f(l::AbstractLattice, styles...)
    for style in styles
        @series l, style
    end
end

@recipe f(l::AbstractLattice, ::Val{:sites}) = @series l, nothing

@recipe function f(l::AbstractLattice, ::Val{:high_contrast})
    markersize := 4
    markercolor := :black
    markerstrokealpha := 1
    markerstrokestyle := :solid
    markerstrokewidth := 2
    markerstrokecolor := :white
    l, nothing
end

@recipe function f(l::AbstractLattice, ::Val{:numbers})
    label --> ""
    annotations = [(" " * string(i), :left, :top, :grey, 6) for i in eachindex(l)]
    @series begin   # The sites
        seriestype := :scatter
        markershape := :none
        series_annotations := annotations
        l, nothing
    end
end

@recipe function f(l::AbstractLattice, ::Val{:bonds}; bonds = (1, 2, 3))
    label := ""
    1 in bonds && @series begin   # The bonds
        seriestype := :path
        label := ""
        l, NearestNeighbor(1)
    end
    2 in bonds && @series begin   # The nnbonds
        seriestype := :path
        linestyle := :dash
        linealpha := 0.5
        label := ""
        l, NearestNeighbor(2)
    end
    3 in bonds && @series begin   # The nnnbonds
        seriestype := :path
        linestyle := :dot
        linealpha := 0.5
        label := ""
        l, NearestNeighbor(3)
    end
end

@recipe function f(l::AbstractLattice, ::Val{:boundaries})
    label := ""
    trs = adapt_boundaries(getboundaries(l), UndefinedLattice())
    seriescolor --> :lightblue
    bonds --> ()
    xmi, xma = extrema(site -> site.coords[1], l)
    xlims --> (xmi + xma) / 2 .+ (xma - xmi) .* (-1, 1)
    ymi, yma = extrema(site -> site.coords[2], l)
    ylims --> (ymi + yma) / 2 .+ (yma - ymi) .* (-1, 1)
    for cind in cartesian_indices(l)
        tup = Tuple(cind)
        all(==(0), tup) && continue
        l2 = l
        for i in eachindex(tup)
            tr = trs[i].translation
            tup[i] < 0 && (tr = -tr)
            for _ in 1:abs(tup[i])
                l2 += tr
            end
        end
        @series begin
            seriesalpha := 0.2
            l2, :sites, :bonds
        end
    end
end

@recipe function f(l::AbstractLattice; showboundaries=true, showbonds=false, shownumbers=false)
    label --> ""
    shownumbers --> showboundaries
    bonds --> (showbonds ? (1,2,3) : ())
    @series l, :bonds
    @series l, :sites
    showboundaries && @series l, :boundaries
    shownumbers && @series l, :numbers
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
    lv.latt, lv.values
end

@recipe function f(ag::AbstractBonds)
    aspect_ratio := :equal
    l = lattice(ag)
    pts = NTuple{dims(l), Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for (s1, s2) in ag
        R = s2.old_site.coords - s1.old_site.coords
        A = s1.site.coords
        B = s2.site.coords
        push!(pts, Tuple(A), Tuple(A + R / 2), br_pt, Tuple(B - R / 2), Tuple(B), br_pt)
    end
    label := nothing
    pts
end

@recipe function f(l::AbstractLattice, b::AbstractBonds)
    adapt_bonds(b, l)
end
