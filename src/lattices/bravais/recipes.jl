using RecipesBase

@recipe function f(l::BravaisLattice, ::Val{:bonds}; bonds = (1, 2, 3))
    label := ""
    1 in bonds && @series begin   # The bonds
        seriestype := :path
        label := ""
        l, default_bonds(l)
    end
    2 in bonds && @series begin   # The nnbonds
        seriestype := :path
        linestyle := :dash
        linealpha := 0.5
        label := ""
        l, default_bonds(l, Val(2))
    end
    3 in bonds && @series begin   # The nnnbonds
        seriestype := :path
        linestyle := :dot
        linealpha := 0.5
        label := ""
        l, default_bonds(l, Val(3))
    end
end

@recipe function f(l::BravaisLattice{N}, ::Val{:boundaries}) where {N}
    label := ""
    for cind in cartesian_indices(l)
        tup = Tuple(cind)
        all(==(0), tup) && continue
        tr_vec = @SVector zeros(Int, N)
        for i in eachindex(tup)
            tr_vec += tup[i] * l.boundaries.bcs[i].R
        end
        shifted_pointers = [shift_site(tr_vec, lp) for lp in l.pointers]
        shifted_lattice = BravaisLattice(l.bravais, shifted_pointers)
        @series begin
            seriesalpha := 0.2
            sites(shifted_lattice)
        end
    end
end

@recipe function f(l::BravaisLattice; showboundaries=false, showbonds=false)
    shownumbers --> showboundaries
    @series sites(l)
    !showbonds && (bonds --> ())
    @series l, :bonds
    showboundaries && @series l, :boundaries
end

function tr_vector(l::BravaisLattice, hop::BravaisShift{<:Pair})
    i, j = hop.site_indices
    return l.bravais.basis[:, j] - l.bravais.basis[:, i] +
     mm_assuming_zeros(l.bravais.translation_vectors, hop.translate_uc)
end
tr_vector(l::BravaisLattice, hop::BravaisShift{Nothing}) =
    mm_assuming_zeros(l.bravais.translation_vectors, hop.translate_uc)
@recipe function f(l::BravaisLattice{N, B}, bss::NTuple{M, BravaisShift} where M) where {N, B}
    aspect_ratio := :equal
    pts = NTuple{N, Float64}[]
    br_pt = fill(NaN, dims(l)) |> Tuple
    for bs in bss
        T = tr_vector(l, bs)
        for i in eachindex(l)
            site1 = l[i]
            p = shift_site(l, site1 + bs)
            p === nothing && continue
            site2 = p[2]
            site2 in l || continue

            a = site1.coords
            b = site2.coords
            push!(pts, Tuple(a), Tuple(a + T / 2), br_pt, Tuple(b), Tuple(b - T / 2), br_pt)
        end
    end
    label := nothing
    pts
end
@recipe f(l::BravaisLattice, bs::BravaisShift, bss::BravaisShift...) = (l, (bs, bss...))

const BravaisLatticeValue{Sym, N} = LatticeValue{<:Number, <:BravaisLattice{N, <:UnitCell{Sym, N}}}
raw"""
    rectified_values(lv::LatticeValue)

Returng an array of the values of `lv` on its macrocell.
The $i$-th element of the array corresponds to the $i$-th site of the macrocell.
If the element is `NaN`, it means that the corresponding site is not present in the `lv`'s lattice.

This function might be quite useful in custom plot recipes.
"""
function rectified_values(lv::BravaisLatticeValue{Sym, N} where Sym) where N
    s = sites(lv)
    mins = Vector(s[1].lp.unit_cell)
    maxs = Vector(s[1].lp.unit_cell)
    for lp in s.latt.pointers
        @. mins = min(mins, lp.unit_cell)
        @. maxs = max(maxs, lp.unit_cell)
    end
    newvals = fill(NaN, Tuple(@. maxs - mins + 1))
    smins = SVector{N}(mins) .- 1
    for (i, lp) in enumerate(s.latt.pointers)
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
    l = sites(lv).latt
    new_l = BravaisLattice(UnitCell{:plot_fallback}(l.bravais.translation_vectors, l.bravais.basis),
        l.pointers)
    LatticeValue(new_l, lv.values)
end

@recipe function f(lv::BravaisLatticeValue{:square, N}) where N
    N < 3 && (seriestype --> :heatmap)
    unitcell = sites(lv).latt.bravais
    if plotattributes[:seriestype] === :heatmap
        axes, heatmap_values = rectified_values(lv)
        c_axes = Tuple(axes[i] .+ unitcell.basis[i, 1] for i in 1:N)
        aspect_ratio := :equal
        c_axes..., transpose(heatmap_values)
    else
        plot_fallback(lv)
    end
end
@recipe function f(lv::LatticeValue{<:Number, <:Sites})
    LatticeValue(lattice(lv).latt, lv.values)
end
