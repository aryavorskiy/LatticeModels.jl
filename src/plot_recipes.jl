@recipe function f(site::LatticeSite)
    seriestype := :scatter
    [Tuple(site.coords),]
end

@recipe function f(l::Lattice; pretty=true, high_contrast=false)
    if high_contrast
        pretty = false
        markersize := 4
        markercolor := :black
        markerstrokealpha := 1
        markerstrokestyle := :solid
        markerstrokewidth := 2
        markerstrokecolor := :white
    end
    if pretty
        l_outp = copy(l)
        fill!(l_outp.mask, true)
        annotations = repeat(Any[""], length(l_outp))
        annotations[l.mask] .= ((i, :left, :top, :grey, 8) for i in 1:length(l))
        series_annotations := annotations
        seriesalpha := l.mask .* 0.9 .+ 0.1
        label --> ""
        l_outp, nothing
    else
        l, nothing
    end
end

@recipe function f(l::Lattice, v)
    aspect_ratio := :equal
    marker_z := v
    pts = collect_coords(l)
    if dims(l) == 3
        X, Y, Z = eachrow(pts)
        Xr, Yr, Zr = eachrow(round.(pts, digits=3))
        seriestype := :scatter3d
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ", ", Zr, ")")
        end
        X, Y, Z
    else
        if dims(l) == 1
            X = vec(pts)
            Y = zero(X)
        else
            X, Y = eachrow(pts[1:2, :])
        end
        seriestype --> :scatter
        if v !== nothing && RecipesBase.is_key_supported(:hover)
            Xr, Yr = eachrow(round.(pts, digits=3))
            hover := string.(round.(v, digits=3), " @ (", Xr, ", ", Yr, ")")
        end
        if plotattributes[:seriestype] == :scatter
            X, Y
        elseif plotattributes[:seriestype] == :surface
            X, Y, v
        else
            throw(ArgumentError("unsupported series type $(plotattributes[:seriestype])"))
        end
    end
end

@recipe function f(lv::PlottableLatticeValue)
    lv.lattice, lv.values
end

@recipe function f(bs::PairSet)
    aspect_ratio := :equal
    l = bs.lattice
    pts = Tuple{Float64,Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for i in 1:length(l)
        site1 = bs.lattice[i]
        A = site1.coords
        for j in 1:length(l)
            if i != j && bs.bmat[i, j]
                site2 = bs.lattice[j]
                B = site2.coords
                T = radius_vector(l, site2, site1)
                push!(pts, Tuple(A))
                push!(pts, Tuple(A + T / 2))
                push!(pts, br_pt)
                push!(pts, Tuple(B))
                push!(pts, Tuple(B - T / 2))
                push!(pts, br_pt)
            end
        end
    end
    label := nothing
    pts
end

@recipe function f(curr::AbstractCurrents)
    l = lattice(curr)
    dims(l) != 2 && error("2D lattice expected")
    Xs = Float64[]
    Ys = Float64[]
    Qs = NTuple{2,Float64}[]
    arrows_scale --> 1
    arrows_rtol --> 1e-2
    seriestype := :quiver
    for (i, site1) in enumerate(l)
        for (j, site2) in enumerate(l)
            j ≥ i && continue
            ij_curr = curr[i, j]::Real
            crd = ij_curr > 0 ? site1.coords : site2.coords
            vc = radius_vector(l, site2, site1)
            vc_n = norm(vc)
            if vc_n < abs(ij_curr * plotattributes[:arrows_scale] / plotattributes[:arrows_rtol])
                push!(Xs, crd[1])
                push!(Ys, crd[2])
                push!(Qs, Tuple(vc * (ij_curr * plotattributes[:arrows_scale] / vc_n)))
            end
        end
    end
    quiver := Qs
    Xs, Ys
end

@recipe function f(tseq::TimeSequence)
    tseq.times, tseq.snapshots
end

@recipe function f(lv::PlottableLatticeValue{:square})
    seriestype --> :heatmap
    if plotattributes[:seriestype] === :heatmap
        aspect_ratio := :equal
        axes_lims = [1:ax for ax in size(lattice(lv))]
        heatmap_values = reshape(macro_cell_values(lv), reverse(size(lattice(lv))))'
        axes_lims..., heatmap_values
    else
        plot_fallback(lv)
    end
end
