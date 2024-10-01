using RecipesBase, ColorTypes

@recipe function f(site::AbstractSite)
    seriestype := :scatter
    [Tuple(site.coords),]
end

@recipe function f(lat::AbstractLattice, style::Symbol)
    lat, Val(style)
end

@recipe function f(::AbstractLattice, ::Val{StyleT}) where StyleT
    error("Unsupported lattice plot style '$StyleT'")
end

function _process_axis(ax)
    if ax isa Int
        return ax
    elseif ax isa Symbol
        p = SitePropertyAlias{ax}()
        if p isa Coord
            return p.axis
        end
    end
    throw(ArgumentError("Invalid axis $ax"))
end
function _get_axes(lat::AbstractLattice, plotattributes)
    lat isa UndefinedLattice && throw(ArgumentError("Cannot plot an undefined lattice."))
    u_axes = plotattributes === nothing ? nothing : get(plotattributes, :axes, nothing)
    if !(u_axes isa Tuple)
        u_axes = (u_axes,)
    end
    length(u_axes) ≤ 3 || throw(ArgumentError("Only ≤3 spatial dimensions are supported."))
    if u_axes !== (nothing,)
        axes = u_axes
    else
        dims(lat) > 3 && throw(ArgumentError("No axes specified for $(dims(lat))D lattice"))
        axes = (:x, :y, :z)[1:dims(lat)]
    end
    axis_numbers = _process_axis.(axes)
    allunique(axis_numbers) || throw(ArgumentError("Duplicate axes in $axes"))
    return axes, axis_numbers
end
@recipe function f(lat::AbstractLattice, ::Val{:sites})
    aspect_ratio := :equal
    seriestype := :scatter
    axes, axis_numbers = _get_axes(lat, plotattributes)
    xguide --> axes[1]
    xwiden --> 1.1
    yguide --> axes[2]
    ywiden --> 1.1
    if length(axes) > 2
        zguide --> axes[3]
        zwiden --> 1.1
    end
    crd = collect_coords(lat)
    @series Tuple(@view(crd[i, :]) for i in axis_numbers)
end

@recipe function f(lat::AbstractLattice, ::Val{:high_contrast})
    label --> :none
    markersize := 4
    markercolor := :black
    markerstrokealpha := 1
    markerstrokestyle := :solid
    markerstrokewidth := 2
    markerstrokecolor := :white
    lat, :sites
end

@recipe function f(lat::AbstractLattice, ::Val{:numbers})
    label --> ""
    inds = get(plotattributes, :site_indices, eachindex(lat))
    alphas = get(plotattributes, :markeralpha, ones(length(lat)))
    annotations = [(" " * string(inds[j]), :left, :top, RGBA(0.5, 0.5, 0.5, alphas[j]), 6)
        for j in eachindex(inds)]
    @series begin   # The sites
        seriestype := :scatter
        markershape := :none
        markeralpha := 0
        series_annotations := annotations
        lat, :sites
    end
end

@recipe function f(lat::AbstractLattice, ::Val{:bonds})
    label := ""
    seriestype := :path
    lss = (:solid, :dash, :dot)
    las = (1, 0.6, 0.5)
    showbonds = get(plotattributes, :showbonds, _showbonds_default(lat))
    if showbonds isa Bool
        showbonds = (showbonds ? _showbonds_default(lat) : ())
    elseif showbonds isa Union{Int, Tuple}
        showbonds = showbonds
    else
        throw(ArgumentError("Invalid `showbonds` argument: " *
            "expected Bool, Int, or NTuple{Int}, got $(typeof(showbonds))"))
    end
    for i in eachindex(showbonds)
        @series begin
            label := ""
            linestype --> lss[min(i, length(lss))]
            linealpha --> las[min(i, length(las))]
            lat, NearestNeighbor(showbonds[i])
        end
    end
end

latticedist2(lat::AbstractLattice, site::AbstractSite) =
    minimum(s -> sum(abs2, s.coords - site.coords), lat)
@recipe f(::AbstractLattice, ::Val{:boundaries}) = @series begin
    aspect_ratio := :equal
    ()
end
@recipe function f(lat::LatticeWithMetadata, ::Val{:boundaries})
    label := ""
    trs = adapt_boundaries(getboundaries(lat), UndefinedLattice())
    for cind in cartesian_indices(lat)
        tup = Tuple(cind)
        all(==(0), tup) && continue
        lat2 = copy(lat)
        for i in eachindex(tup)
            tr = trs[i].translation
            tup[i] < 0 && (tr = -tr)
            for _ in 1:abs(tup[i])
                lat2 += tr
            end
        end
        ldists = [latticedist2(lat, site) for site in lat2]
        diam2 = maximum(ldists, init=1.0)
        q = max(1/3.1, 3.5/√diam2)
        is = findall(<(diam2 * q^2), ldists)
        site_indices := is
        alphafalloff = map(x -> 0.65 * exp(-2.5x / (diam2 * q^2)), @view ldists[is])
        @series begin
            sitealpha := alphafalloff
            seriescolor --> :lightblue
            lat2[is], :bonds
        end
        @series begin
            markeralpha := alphafalloff
            markerstrokewidth := 0.5
            if :marker_z in keys(plotattributes)
                marker_z := plotattributes[:marker_z][is]
            else
                seriescolor --> :lightblue
            end
            stripmeta(lat2)[is]
        end
    end
end

_showbonds_default(::AbstractLattice) = ()
_showbonds_default(l::LatticeWithMetadata) =
    length(getnnbonds(l)) ≥ 1 ? 1 : _showbonds_default(stripmeta(l))
@recipe function f(lat::AbstractLattice; showboundaries=true, shownumbers=false)
    label --> ""
    @series lat, :bonds
    @series lat, :sites
    showboundaries && @series lat, :boundaries
    shownumbers && @series lat, :numbers
end

function getcircle(r::Real, n::Int)
    is = (0:n) .* 2pi/n
    return r * cos.(is), r * sin.(is)
end
function getshape(lat::AbstractLattice)
    site = lat[1]
    r = sitedistance(lat, site, lat[2])
    for i in 3:length(lat)
        r = min(r, sitedistance(lat, site, lat[i]))
    end
    return getcircle(r/2, 20)
end
getshape(lat::AbstractLattice, _) = getshape(stripmeta(lat))   # fallback for undefined lattice types
function getshape(lat::LatticeWithMetadata)
    if hasmeta(lat, :latticetype)
        return getshape(stripmeta(lat), gettype(lat))
    elseif length(getnnbonds(lat)) > 0
        nnb = getnnbonds(lat)
        return getcircle(nnb.dists[1], 20)
    else
        return getshape(stripmeta(lat))
    end
end
function moveshape(shape, loc::SVector{2}, scale=1)
    xs, ys = shape
    return xs .* scale .+ loc[1], ys .* scale .+ loc[2]
end
@recipe function f(lv::LatticeValue{<:Number}, ::Val{:tiles})
    dims(lv) == 2 || throw(ArgumentError("Only 2D lattices are supported for tiles plot"))
    :axes in keys(plotattributes) && throw(ArgumentError("Cannot use `axes` argument with tiles plot"))
    scale_markers = get(plotattributes, :markerscale, false)
    @assert scale_markers isa Real "`markerscale` must be a real number or `Bool`"
    mx = maximum(abs, lv.values)
    lat = lattice(lv)

    seriescolor --> :matter
    shapetype = get(plotattributes, :markershape, :polygon)
    if shapetype === :circle
        shape = getshape(lat, nothing)  # Trigger fallback
    elseif shapetype === :polygon
        shape = getshape(lat)
    else throw(ArgumentError("Unsupported shape type `:$shapetype`"))
    end
    if scale_markers !== false
        @series begin
            seriestype := :path
            markershape := :none
            linecolor := :grey
            linealpha := 0.5
            linewidth := 2
            lat, :bonds
        end
    end
    xguide --> "x"
    yguide --> "y"
    for site in lv.lat
        @series begin
            seriestype := :shape
            markershape := :none
            fill_z := lv[site]
            aspect_ratio := :equal
            moveshape(shape, site.coords,
                scale_markers === false ? 1 : lv[site] / mx * scale_markers)
        end
    end
end
@recipe function f(lv::LatticeValue{<:Number}, ::Val{:scatter})
    scale_markers = get(plotattributes, :markerscale, true)
    @assert scale_markers isa Real "`markerscale` must be a real number or `Bool`"
    mx = maximum(abs, lv.values)
    lat = lattice(lv)

    seriescolor --> :matter
    linecolor := :grey
    linealpha := 0.5
    linewidth := 2
    marker_z := lv.values
    marker_sz = get(plotattributes, :markersize, scale_markers === false ? 4 : 8)
    @series begin
        markersize := 0
        (lat,)
    end
    @series begin
        markersize := scale_markers === false ? marker_sz :
            @. marker_sz * abs(lv.values) / mx * scale_markers
        markerstrokewidth --> 0.5
        aspect_ratio := :equal
        lat, :sites
    end
end
@recipe function f(lv::LatticeValue{<:Number}, ::Val{:line})
    seriestype --> :path
    axis, i = only.(_get_axes(lattice(lv), plotattributes))
    l = lattice(lv)
    i in 1:dims(l) || throw(ArgumentError("Invalid axis $axis"))
    @series begin
        xguide --> i in 1:3 ? (:x, :y, :z)[i] : raw"$x_{$i}"
        crd = collect_coords(l)[i, :]
        perm = sortperm(crd)
        crd[perm], lv.values[perm]
    end
end

function heatmap_data(lv::LatticeValue{T}, axis_numbers, bins) where {T<:Real}
    crd = collect_coords(lattice(lv))[collect(axis_numbers), :]
    min_pt = dropdims(minimum(crd, dims=2), dims=2)
    max_pt = dropdims(maximum(crd, dims=2), dims=2)
    if bins isa Number
        xbins, ybins = round.(Int, sqrt(bins / prod(max_pt - min_pt)) * (max_pt - min_pt))
    elseif bins isa NTuple{2, Any}
        xbins, ybins = bins
    else
        throw(ArgumentError("Invalid `bins` argument: expected integer or 2-tuple, got $(typeof(bins))"))
    end
    if xbins === ybins === nothing
        r = sqrt(prod(max_pt - min_pt))
        xbins, ybins = round.(Int, (max_pt - min_pt) / r * sqrt(length(lv)) / √2)
    elseif xbins === nothing
        xbins = round(Int, ybins * (max_pt[1] - min_pt[1]) / (max_pt[2] - min_pt[2]))
    elseif ybins === nothing
        ybins = round(Int, xbins * (max_pt[2] - min_pt[2]) / (max_pt[1] - min_pt[1]))
    end
    bins = (xbins, ybins)
    prod(bins) > length(lv) &&
        @warn "$xbins×$ybins ($(prod(bins))) is too many bins for $(length(lv)) data points"
    counts = zeros(Int, bins)
    sums = zeros(float(T), bins)
    @inbounds @simd for i in Base.axes(crd, 2)
        ix = round(Int, (crd[1, i] - min_pt[1]) / (max_pt[1] - min_pt[1]) * (xbins - 1)) + 1
        iy = round(Int, (crd[2, i] - min_pt[2]) / (max_pt[2] - min_pt[2]) * (ybins - 1)) + 1
        counts[ix, iy] += 1
        sums[ix, iy] += lv.values[i]
    end
    xs = range(min_pt[1], max_pt[1], length=xbins)
    ys = range(min_pt[2], max_pt[2], length=ybins)
    return xs, ys, transpose(sums ./ counts)
end

@recipe function f(lv::LatticeValue{T}, ::Val{:hmap}; xbins=nothing, ybins=nothing, bins=(xbins, ybins)) where {T<:Number}
    aspect_ratio := :equal
    axes, axis_numbers = _get_axes(lattice(lv), plotattributes)
    if length(axes) != 2
        throw(ArgumentError("2D axes expected; got $axes"))
    end
    xguide --> axes[1]
    yguide --> axes[2]
    if plotattributes[:seriestype] == :histogram2d
        plotattributes[:seriestype] = :heatmap
    elseif plotattributes[:seriestype] == :contour
        fill --> true
        linewidth --> 0.7
    else
        error("Unsupported seriestype $(plotattributes[:seriestype]) for hmap plot")
    end
    heatmap_data(lv, axis_numbers, bins)
end

@recipe function f(lv::LatticeValue{<:Number})
    label --> ""
    axes, _ = _get_axes(lattice(lv), plotattributes)
    if length(axes) == 1
        @series lv, Val(:line)
    else
        seriestype --> :scatter
        if plotattributes[:seriestype] in (:shape, :heatmap) && dims(lv) == 2
            @series lv, Val(:tiles)
        elseif plotattributes[:seriestype] in (:histogram2d, :contour)
            @series lv, Val(:hmap)
        elseif plotattributes[:seriestype] == :scatter
            @series lv, Val(:scatter)
        else
            error("Unsupported seriestype $(plotattributes[:seriestype])")
        end
    end
end

function collect_coords(bonds::AbstractBonds, sitealpha=nothing)
    lat = lattice(bonds)
    pts = SVector{dims(lat), Float64}[]
    alphas = Float64[]
    nans = zero(SVector{dims(lat)}) * NaN
    for (s1, s2) in bonds
        A = s1.site.coords
        B = s2.site.coords
        R = s2.old_site.coords - s1.old_site.coords
        if R ≈ B - A
            push!(pts, A, B, nans)
            if sitealpha !== nothing
                push!(alphas, sitealpha[s1.index], sitealpha[s1.index], NaN)
            end
        else
            push!(pts, A, A + R, nans, B - R, B, nans)
            if sitealpha !== nothing
                push!(alphas, sitealpha[s1.index], sitealpha[s1.index], NaN,
                    sitealpha[s2.index], sitealpha[s2.index], NaN)
            end
        end
    end
    zs = zeros(dims(lat), length(pts))
    for i in 1:length(pts)
        zs[:, i] = pts[i]
    end
    return zs, alphas
end

@recipe function f(bonds::AbstractBonds)
    aspect_ratio := :equal
    label := nothing
    axes, axis_numbers = _get_axes(lattice(bonds), plotattributes)
    xguide --> axes[1]
    yguide --> axes[2]
    if length(axes) > 2
        zguide --> axes[3]
    end
    crd, linealphas = collect_coords(bonds, get(plotattributes, :sitealpha, nothing))
    if !isempty(linealphas)
        linealpha := linealphas[1:end-1]
    end
    @series Tuple(@view(crd[i, :]) for i in axis_numbers)
end

@recipe function f(l::AbstractLattice, b::AbstractBonds)
    adapt_bonds(b, l)
end
