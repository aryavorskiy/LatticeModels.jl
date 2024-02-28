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
function _get_axes(l, default)
    l isa UndefinedLattice && throw(ArgumentError("Cannot plot an undefined lattice."))
    axes = default !== nothing ? default :
        (dims(l) == 3 ? (:x, :y, :z) : (:x, :y))
    if length(axes) < 2
        throw(ArgumentError("At least two axes are required; got $axes"))
    elseif length(axes) > dims(l)
        throw(ArgumentError("Cannot display a $(dims(l))D lattice in $(length(axes)) axes"))
    end
    axis_numbers = _process_axis.(axes)
    allunique(axis_numbers) || throw(ArgumentError("Duplicate axes in $axes"))
    return axes, axis_numbers
end
@recipe function f(l::AbstractLattice, ::Val{:sites})
    aspect_ratio := :equal
    seriestype := :scatter
    axes, axis_numbers = _get_axes(l, get(plotattributes, :axes, nothing))
    xguide --> axes[1]
    yguide --> axes[2]
    if length(axes) > 2
        zguide --> axes[3]
    end
    rows = eachrow(collect_coords(l))
    @series Tuple(rows[i] for i in axis_numbers)
end

@recipe function f(l::AbstractLattice, ::Val{:high_contrast})
    markersize := 4
    markercolor := :black
    markerstrokealpha := 1
    markerstrokestyle := :solid
    markerstrokewidth := 2
    markerstrokecolor := :white
    l, :sites
end

@recipe function f(l::AbstractLattice, ::Val{:numbers})
    label --> ""
    annotations = [(" " * string(i), :left, :top, :grey, 6) for i in eachindex(l)]
    @series begin   # The sites
        seriestype := :scatter
        markershape := :none
        markeralpha := 0
        series_annotations := annotations
        l, :sites
    end
end

@recipe function f(l::AbstractLattice, ::Val{:bonds})
    label := ""
    seriestype := :path
    lss = (:solid, :dash, :dot)
    las = (1, 0.6, 0.5)
    showbonds = get(plotattributes, :showbonds, 1)
    for i in eachindex(showbonds)
        @series begin
            label := ""
            linestype := lss[min(i, length(lss))]
            linealpha := las[min(i, length(las))]
            l, NearestNeighbor(showbonds[i])
        end
    end
end

@recipe function f(l::AbstractLattice, ::Val{:boundaries})
    label := ""
    trs = adapt_boundaries(getboundaries(l), UndefinedLattice())
    seriescolor --> :lightblue
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

@recipe function f(l::AbstractLattice; showboundaries=true, showbonds=true, shownumbers=false)
    label --> ""
    shownumbers --> showboundaries
    if showbonds isa Bool
        showbonds := (showbonds ? (1,) : ())
    elseif showbonds isa Union{Int, Tuple}
        showbonds := showbonds
    else
        throw(ArgumentError(""))
    end
    @series l, :bonds
    @series l, :sites
    showboundaries && @series l, :boundaries
    shownumbers && @series l, :numbers
end

function getcircle(r::Real, n::Int)
    is = (0:n) .* 2pi/n
    return r * cos.(is), r * sin.(is)
end
function getshape(l::AbstractLattice)
    site = l[1]
    r = site_distance(l, site, l[2])
    for i in 3:length(l)
        r = min(r, site_distance(l, site, l[i]))
    end
    return getcircle(r/2, 20)
end
getshape(l::AbstractLattice, _) = getshape(l)   # fallback for undefined lattice types
function getshape(l::LatticeWithParams)
    if hasparam(l, :latticetype)
        return getshape(l.lat, l.latticetype)
    elseif length(getnnbonds(l)) > 0
        nnb = getnnbonds(l)
        return getcircle(nnb.dists[1], 20)
    else
        return getshape(l.lat)
    end
end
function moveshape(shape, loc::SVector{2})
    xs, ys = shape
    return xs .+ loc[1], ys .+ loc[2]
end
@recipe function f(lv::LatticeValue{<:Number})
    seriestype --> :image
    if plotattributes[:seriestype] !== :image || dims(lattice(lv)) != 2
        marker_z := lv.values
        lv.latt, :sites
    else
        label := ""
        aspect_ratio := :equal
        shape = getshape(lv.latt)
        for site in lv.latt
            @series begin
                seriestype := :shape
                fill_z := lv[site]
                moveshape(shape, site.coords)
            end
        end
    end
end

function collect_coords(bonds::AbstractBonds)
    l = lattice(bonds)
    pts = SVector{dims(l), Float64}[]
    nans = fill(NaN, dims(l)) |> Tuple |> SVector
    for (s1, s2) in bonds
        A = s1.site.coords
        B = s2.site.coords
        R = s2.old_site.coords - s1.old_site.coords
        if R ≈ B - A
            push!(pts, A, B, nans)
        else
            push!(pts, A, A + R / 2, nans, B - R / 2, B, nans)
        end
    end
    zs = zeros(dims(l), length(pts))
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
