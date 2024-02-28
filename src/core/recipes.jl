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
    seriescolor --> :lightblue
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
        alphafalloff = map(x -> 0.7 * exp(-2.5x / (diam2 * q^2)), @view ldists[is])
        @series begin
            markeralpha := alphafalloff
            lat2[is], :sites
        end
        @series begin
            linealpha := 0.5
            lat2[is], :bonds
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
    r = site_distance(lat, site, lat[2])
    for i in 3:length(lat)
        r = min(r, site_distance(lat, site, lat[i]))
    end
    return getcircle(r/2, 20)
end
getshape(lat::AbstractLattice, _) = getshape(lat)   # fallback for undefined lattice types
function getshape(lat::LatticeWithParams)
    if hasparam(lat, :latticetype)
        return getshape(stripparams(lat), lat.latticetype)
    elseif length(getnnbonds(lat)) > 0
        nnb = getnnbonds(lat)
        return getcircle(nnb.dists[1], 20)
    else
        return getshape(stripparams(lat))
    end
end
function moveshape(shape, loc::SVector{2})
    xs, ys = shape
    return xs .+ loc[1], ys .+ loc[2]
end
@recipe function f(lv::LatticeValue{<:Number}; showbonds::Bool=true)
    seriestype --> (length(lv) < 500  ? :image : :scatter)
    seriescolor --> :tempo
    label --> ""
    lat = lattice(lv)
    if plotattributes[:seriestype] === :image && dims(lv) == 2
        label := ""
        aspect_ratio := :equal
        shape = getshape(lat)
        for site in lv.latt
            @series begin
                seriestype := :shape
                fill_z := lv[site]
                moveshape(shape, site.coords)
            end
        end
    else
        showbonds && @series begin
            showbonds := 1
            seriescolor := :grey
            linealpha := 0.5
            linewidth := 2
            lat, :bonds
        end
        marker_z := lv.values
        markersize := 5
        @series lat, :sites
    end
end

function collect_coords(bonds::AbstractBonds)
    lat = lattice(bonds)
    pts = SVector{dims(lat), Float64}[]
    nans = fill(NaN, dims(lat)) |> Tuple |> SVector
    for (s1, s2) in bonds
        A = s1.site.coords
        B = s2.site.coords
        R = s2.old_site.coords - s1.old_site.coords
        if R ≈ B - A
            push!(pts, A, B, nans)
        else
            push!(pts, A, A + R, nans, B - R, B, nans)
        end
    end
    zs = zeros(dims(lat), length(pts))
    for i in 1:length(pts)
        zs[:, i] = pts[i]
    end
    return zs
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
    rows = eachrow(collect_coords(bonds))
    @series Tuple(rows[i] for i in axis_numbers)
end

@recipe function f(l::AbstractLattice, b::AbstractBonds)
    adapt_bonds(b, l)
end
