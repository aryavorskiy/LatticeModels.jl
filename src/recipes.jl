function dummy_lattice(uc::UnitCell{N,NU,NB}) where {N,NU,NB}
    l = BravaisLattice(uc, [BravaisPointer(zero(SVector{NU,Int}), i) for i in 1:NB])
    for j in 1:NU
        vc = SVector{NU,Int}(eachindex(zero(SVector{NU,Int})) .== j)
        for k in 1:NB
            push!(l, BravaisPointer(vc, k), BravaisPointer(-vc, k))
        end
    end
    return l
end

@recipe function f(uc::UnitCell{N,NU,NB}; offset=:origin, rotate=nothing) where {N,NU,NB}
    uc = transform_unitcell(uc; offset=offset, rotate=rotate)
    l = dummy_lattice(uc)
    markeralphas_one = fill(0.5, 1 + 2NU)
    markeralphas_one[NU + 1] = 1
    @series begin
        linecolor := :grey
        linealpha := 0.5
        l, :bonds
    end
    for i in 1:NB
        @series begin
            seriestype := :scatter
            markersize --> 15
            markeralpha --> markeralphas_one
            aspect_ratio := :equal
            ywiden := 1.3
            label := ""
            l[BasisIndex() => i], :sites
        end
    end
    @series begin
        seriestype := :quiver
        linewidth := 2
        aspect_ratio := :equal
        ywiden := 1.2
        quiver := [Tuple(unitvector(uc, j)) for j in 1:NU]
        zeros(NU), zeros(NU)
    end
end

@recipe function f(tseq::TimeSequence)
    xlabel --> "t"
    label --> ""
    tseq.times, tseq.values
end

@recipe function f(curr::AbstractCurrents;
        arrowheadsize=0.15, arrowheadwidth=1/3, arrowtransparency=true)
    lat = lattice(curr)
    axes, axis_numbers = _get_axes(lat, get(plotattributes, :axes, nothing))
    length(axes) != 2 && error("2D axes expected; got $(axes)D")
    ns = SVector(axis_numbers)
    xguide --> axes[1]
    yguide --> axes[2]
    label --> :none
    Xs = Float64[]
    Ys = Float64[]
    @inline function _pushpts!(pts...)
        for pt in pts
            push!(Xs, pt[1])
            push!(Ys, pt[2])
        end
        push!(Xs, NaN)
        push!(Ys, NaN)
    end
    Zs = Float64[]
    Is, Js, Vs = findnz(curr)
    for ind in 1:length(Is)
        site1 = lat[Is[ind]]
        site2 = lat[Js[ind]]
        val = Vs[ind]
        if val < 0
            site1, site2 = site2, site1
            val = -val
        end
        v1 = site1.coords[ns]
        v2 = site2.coords[ns]
        d = normalize(v2 - v1)
        o = SVector(d[2], -d[1])
        _pushpts!(v1, v2, v2 - arrowheadsize * (d - o * arrowheadwidth),
            v2 - arrowheadsize * (d + o * arrowheadwidth), v2)
        push!(Zs, val, val, val, val, val, NaN)
    end
    if isempty(Zs)
        _pushpts!()
        push!(Zs, NaN)
    end
    @series begin
        seriestype := :path
        seriescolor --> :matter
        linewidth --> 2.5
        if arrowtransparency
            mx = maximum(x -> isnan(x) ? zero(x) : x, Zs)
            linealpha --> @view(Zs[1:end-1]) / mx
        end
        line_z --> Zs
        Xs, Ys
    end
    @series begin
        seriestype := :scatter
        markersize := 0.5
        markercolor := :grey
        markeralpha := 0.2
        lattice(curr), :sites
    end
end

@recipe function f(shapes::AbstractShape...;)
    c = get(plotattributes, :scale, nothing)
    seriescolor --> :grey
    linestyle --> :dash
    label --> ""
    for shape in shapes
        if c === nothing
            @series topath2d(shape)
        else
            @series topath2d(scale(shape, c))
        end
    end
end
