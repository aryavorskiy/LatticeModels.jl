using RecipesBase

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

@recipe function f(lat::AbstractLattice, styles...)
    for style in styles
        @series lat, style
    end
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
function _get_axes(lat, default)
    lat isa UndefinedLattice && throw(ArgumentError("Cannot plot an undefined lattice."))
    axes = default !== nothing ? default :
        (dims(lat) == 3 ? (:x, :y, :z) : (:x, :y))
    if length(axes) < 2
        throw(ArgumentError("At least two axes are required; got $axes"))
    elseif length(axes) > dims(lat)
        throw(ArgumentError("Cannot display a $(dims(lat))D lattice in $(length(axes)) axes"))
    end
    axis_numbers = _process_axis.(axes)
    allunique(axis_numbers) || throw(ArgumentError("Duplicate axes in $axes"))
    return axes, axis_numbers
end
@recipe function f(lat::AbstractLattice, ::Val{:sites})
    aspect_ratio := :equal
    seriestype := :scatter
    axes, axis_numbers = _get_axes(lat, get(plotattributes, :axes, nothing))
    xguide --> axes[1]
    yguide --> axes[2]
    if length(axes) > 2
        zguide --> axes[3]
    end
    rows = eachrow(collect_coords(lat))
    @series Tuple(rows[i] for i in axis_numbers)
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
    annotations = [(" " * string(i), :left, :top, :grey, 6) for i in eachindex(lat)]
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
    showbonds = get(plotattributes, :showbonds, 1)
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
@recipe function f(lat::AbstractLattice, ::Val{:boundaries})
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
        diam2 = maximum(ldists)
        q = max(1/3.1, 3.5/√diam2)
        is = findall(<(diam2 * q^2), ldists)
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
            lat2[is], :sites
        end
    end
end

@recipe function f(lat::AbstractLattice; showboundaries=true, showbonds=true, shownumbers=false)
    label --> ""
    shownumbers --> showboundaries
    if showbonds isa Bool
        showbonds := (showbonds ? (1,) : ())
    elseif showbonds isa Union{Int, Tuple}
        showbonds := showbonds
    else
        throw(ArgumentError(""))
    end
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
getshape(lat::AbstractLattice, _) = getshape(stripparams(lat))   # fallback for undefined lattice types
function getshape(lat::LatticeWithParams)
    if hasparam(lat, :latticetype)
        return getshape(stripparams(lat), gettype(lat))
    elseif length(getnnbonds(lat)) > 0
        nnb = getnnbonds(lat)
        return getcircle(nnb.dists[1], 20)
    else
        return getshape(stripparams(lat))
    end
end
function moveshape(shape, loc::SVector{2}, scale=1)
    xs, ys = shape
    return xs .* scale .+ loc[1], ys .* scale .+ loc[2]
end
@recipe function f(lv::LatticeValue{<:Number};)
    seriestype --> :scatter
    seriescolor --> :matter
    label --> ""
    aspect_ratio := :equal
    markerscale --> plotattributes[:seriestype] == :scatter
    scale_markers = plotattributes[:markerscale]
    @assert scale_markers isa Real "`markerscale` must be a real number or `Bool`"
    mx = maximum(abs, lv.values)
    lat = lattice(lv)
    if plotattributes[:seriestype] === :shape && dims(lv) == 2
        :axes in keys(plotattributes) && throw(ArgumentError("Cannot use `:axes` with `:shape` seriestype"))
        shapetype = get(plotattributes, :markershape, :polygon)
        if shapetype === :circle
            shape = getshape(lat, nothing)  # Trigger fallback
        elseif shapetype === :polygon
            shape = getshape(lat)
        else throw(ArgumentError("Unsupported shape type `:$shapetype`"))
        end
        if scale_markers !== :none
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
                moveshape(shape, site.coords,
                    scale_markers === false ? 1 : lv[site] / mx * scale_markers)
            end
        end
    else
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
            markerstrokewidth := 0.5
            lat, :sites
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
    axes, axis_numbers = _get_axes(lattice(bonds), get(plotattributes, :axes, nothing))
    xguide --> axes[1]
    yguide --> axes[2]
    if length(axes) > 2
        zguide --> axes[3]
    end
    crd, linealphas = collect_coords(bonds, get(plotattributes, :sitealpha, nothing))
    rows = eachrow(crd)
    if !isempty(linealphas)
        linealpha := linealphas[1:end-1]
    end
    @series Tuple(rows[i] for i in axis_numbers)
end

@recipe function f(l::AbstractLattice, b::AbstractBonds)
    adapt_bonds(b, l)
end
