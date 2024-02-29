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
    l = dummy_lattice(transform_unitcell(uc; offset=offset, rotate=rotate))
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
            aspect_ratio := 1
            ywiden := 1.3
            label := ""
            l[BasisIndex() => i], :sites
        end
    end
    @series begin
        seriestype := :quiver
        linewidth := 2
        aspect_ratio := 1
        ywiden := 1.2
        quiver := [Tuple(unitvector(uc, j)) for j in 1:NU]
        zeros(NU), zeros(NU)
    end
end

@recipe function f(tseq::TimeSequence)
    tseq.times, tseq.values
end

@recipe function f(curr::AbstractCurrents; showsites=false)
    lat = lattice(curr)
    dims(lat) != 2 && error("2D lattice expected")
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
    Vs = Float64[]
    for ((site1, site2), val) in curr
        abs(val) < 1e-10 && continue
        if val < 0
            site1, site2 = site2, site1
            val = -val
        end
        v1 = site1.coords
        v2 = site2.coords
        d = normalize(v2 - v1)
        o = SVector(d[2], -d[1])
        _pushpts!(v1, v2, v2 - 0.15d - 0.05o, v2 - 0.15d + 0.05o, v2)
        push!(Vs, val, val, val, val, val, NaN)
    end
    @series begin
        seriestype := :path
        seriescolor --> :tempo
        linewidth --> 2.5
        line_z --> Vs
        Xs, Ys
    end
    showsites && @series begin
        seriestype := :scatter
        markersize := 2.0
        markercolor := :gray
        markeralpha := 0.6
        markerstrokewidth := 0
        label := ""
        lattice(curr), :sites
    end
end
