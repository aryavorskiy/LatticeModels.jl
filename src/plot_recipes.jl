using RecipesBase

@recipe function f(site::LatticeSite)
    seriestype := :scatter
    [Tuple(site.coords),]
end

@recipe function f(l::Lattice, style::Symbol=:plain)
    l, Val(style)
end

@recipe function f(::Lattice, ::Val{StyleT}) where StyleT
    error("Unsupported lattice plot style $StyleT")
end

@recipe function f(l::Lattice, ::Val{:plain})
    l, nothing
end

@recipe function f(l::Lattice, ::Val{:high_contrast})
    markersize := 4
    markercolor := :black
    markerstrokealpha := 1
    markerstrokestyle := :solid
    markerstrokewidth := 2
    markerstrokecolor := :white
    l, nothing
end

@recipe function f(l::Lattice, ::Val{:pretty})
    label --> ""
    @series begin   # The sites
        seriestype := :scatter
        annotations = [(" " * string(i), :left, :top, :grey, 6) for i in 1:length(l)]
        series_annotations := annotations
        l, nothing
    end

    @series begin   # The bonds
        seriestype := :path
        label := ""
        l, default_bonds(l)
    end
    @series begin   # The nnbonds
        seriestype := :path
        linestyle := :dash
        linealpha := 0.8
        label := ""
        l, default_bonds(l, Val(2))
    end
    @series begin   # The nnnbonds
        seriestype := :path
        linestyle := :dot
        linealpha := 0.5
        label := ""
        l, default_bonds(l, Val(3))
    end
end

@recipe function f(l::Lattice, v)
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

function displace_site(l::Lattice, site::LatticeSite, bs::SiteOffset)
    new_site = get_site_periodic(l, site + bs)
    new_site in l ? new_site : nothing
end

@recipe function f(ag::AbstractGraph)
    aspect_ratio := :equal
    l = lattice(ag)
    pts = NTuple{dims(l), Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for i in 1:length(l), j in 1:length(l)
        site1 = l[i]
        site2 = l[j]
        !match(ag, site1, site2) && continue
        A = site1.coords
        B = site2.coords
        T = radius_vector(l, site1, site2)
        push!(pts, Tuple(A), Tuple(A + T / 2), br_pt, Tuple(B), Tuple(B - T / 2), br_pt)
    end
    label := nothing
    pts
end

@recipe function f(l::Lattice{N, B}, bss::NTuple{M, SiteOffset} where M) where {N, B}
    aspect_ratio := :equal
    pts = NTuple{N, Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for bs in bss
        T = radius_vector(l, bs)
        for i in 1:length(l)
            site1 = l[i]
            site2 = displace_site(l, site1, bs)
            site2 === nothing && continue

            a = site1.coords
            b = site2.coords
            push!(pts, Tuple(a), Tuple(a + T / 2), br_pt, Tuple(b), Tuple(b - T / 2), br_pt)
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
    tseq.times, tseq.values
end

raw"""
    rectified_values(lv::LatticeValue)

Returng an array of the values of `lv` on its macrocell.
The $i$-th element of the array corresponds to the $i$-th site of the macrocell.
If the element is `NaN`, it means that the corresponding site is not present in the `lv`'s lattice.

This function might be quite useful in custom plot recipes.
"""
function macro_cell_values(lv::LatticeValue{<:Number, <:Lattice{<:Bravais{Sym,N,1}, N} where Sym}) where N
    l = lattice(lv)
    mins = Vector(l[1].unit_cell)
    maxs = Vector(l[1].unit_cell)
    for lp in l.pointers
        @. mins = min(mins, lp.unit_cell)
        @. maxs = max(maxs, lp.unit_cell)
    end
    newvals = fill(NaN, Tuple(@. maxs - mins + 1))
    smins = SVector{N}(mins) .- 1
    for (i, lp) in enumerate(l.pointers)
        newvals[(lp.unit_cell - smins)...] = lv.values[i]
    end
    newvals
end

"""
    plot_fallback(lv::LatticeValue)

Creates a copy of `lv` lattice value with its `LatticeSym` overwritten to `:plot_fallback`.
Use it to invoke the default plot recipe for `LatticeValues` when defining a custom one.
"""
function plot_fallback(lv::LatticeValue)
    l = lattice(lv)
    new_l = Lattice(Bravais(:plot_fallback)(l.bravais.translation_vectors, l.bravais.basis),
        l.pointers)
    LatticeValue(new_l, lv.values)
end

const PlottableLatticeValue{Sym} = LatticeValue{<:Number, <:Lattice{N, Bravais{Sym, N}} where N}

@recipe function f(lv::PlottableLatticeValue{:square})
    seriestype --> :heatmap
    if plotattributes[:seriestype] === :heatmap
        aspect_ratio := :equal
        axes_lims = [1:ax for ax in macrocell_size(lattice(lv))]
        heatmap_values = reshape(macro_cell_values(lv), reverse(macrocell_size(lattice(lv))))'
        axes_lims..., heatmap_values
    else
        plot_fallback(lv)
    end
end
