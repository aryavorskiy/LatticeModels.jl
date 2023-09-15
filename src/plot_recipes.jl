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
        annotations = [(" " * string(i), :left, :top, :grey, 6) for i in 1:length(l)]
        series_annotations := annotations
        l, nothing
    end
end

@recipe function f(l::BravaisLattice, ::Val{:bonds})
    label := ""
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

@recipe function f(l::BravaisLattice, ::Val{:pretty})
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

function displace_site(l::BravaisLattice, site::BravaisSite, bs::SiteOffset)
    new_site = get_site_periodic(l, site + bs)
    new_site in l ? new_site : nothing
end

@recipe function f(ag::AdjacencyMatrix)
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
        push!(pts, Tuple(A), Tuple(B), br_pt)
    end
    label := nothing
    pts
end

function tr_vector(l::BravaisLattice, hop::SiteOffset{<:Pair})
    i, j = hop.site_indices
    return bravais(l).basis[:, j] - bravais(l).basis[:, i] +
     mm_assuming_zeros(bravais(l).translation_vectors, hop.translate_uc)
end
tr_vector(l::BravaisLattice, hop::SiteOffset{Nothing}) =
    mm_assuming_zeros(l.bravais.translation_vectors, hop.translate_uc)
@recipe function f(l::BravaisLattice{N, B}, bss::NTuple{M, SiteOffset} where M) where {N, B}
    aspect_ratio := :equal
    pts = NTuple{N, Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for bs in bss
        T = tr_vector(l, bs)
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
            vc = site2.coords - site1.coords
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

const BravaisLatticeValue{Sym, N} = LatticeValue{<:Number, <:BravaisLattice{N, <:UnitCell{Sym, N}}}
raw"""
    rectified_values(lv::LatticeValue)

Returng an array of the values of `lv` on its macrocell.
The $i$-th element of the array corresponds to the $i$-th site of the macrocell.
If the element is `NaN`, it means that the corresponding site is not present in the `lv`'s lattice.

This function might be quite useful in custom plot recipes.
"""
function rectified_values(lv::BravaisLatticeValue{Sym, N} where Sym) where N
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
    Tuple(mins[i]:maxs[i] for i in 1:N), newvals
end

"""
    plot_fallback(lv::BravaisLatticeValue)

Creates a copy of `lv` lattice value with the underlying `LatticeSym` overwritten to `:plot_fallback`.
Use it to invoke the default plot recipe for `LatticeValues` when defining a custom one.
"""
function plot_fallback(lv::BravaisLatticeValue)
    l = lattice(lv)
    new_l = BravaisLattice(UnitCell(:plot_fallback)(l.bravais.translation_vectors, l.bravais.basis),
        l.pointers)
    LatticeValue(new_l, lv.values)
end

@recipe function f(lv::BravaisLatticeValue{:square, N}) where N
    N < 3 && (seriestype --> :heatmap)
    if plotattributes[:seriestype] === :heatmap
        axes, heatmap_values = rectified_values(lv)
        c_axes = Tuple(axes[i] .+ lattice(lv).bravais.basis[i, 1] for i in 1:N)
        aspect_ratio := :equal
        c_axes..., transpose(heatmap_values)
    else
        plot_fallback(lv)
    end
end
