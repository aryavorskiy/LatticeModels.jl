using RecipesBase

@recipe function f(l::BravaisLattice, ::Val{:bonds}; bonds = (1, 2, 3))
    label := ""
    1 in bonds && @series begin   # The bonds
        seriestype := :path
        label := ""
        l, NearestNeighbor(1)
    end
    2 in bonds && @series begin   # The nnbonds
        seriestype := :path
        linestyle := :dash
        linealpha := 0.5
        label := ""
        l, NearestNeighbor(2)
    end
    3 in bonds && @series begin   # The nnnbonds
        seriestype := :path
        linestyle := :dot
        linealpha := 0.5
        label := ""
        l, NearestNeighbor(3)
    end
end

@recipe function f(l::BravaisLattice{N}, ::Val{:boundaries}) where {N}
    label := ""
    for cind in cartesian_indices(l)
        tup = Tuple(cind)
        all(==(0), tup) && continue
        tr_vec = @SVector zeros(Int, N)
        for i in eachindex(tup)
            tr_vec += tup[i] * getboundaries(l)[i].R
        end
        shifted_pointers = [shift_site(tr_vec, lp) for lp in l.pointers]
        shifted_lattice = BravaisLattice(l.unitcell, shifted_pointers)
        @series begin
            seriesalpha := 0.2
            shifted_lattice, :sites
        end
    end
end

@recipe function f(l::BravaisLattice; showboundaries=false, showbonds=false)
    shownumbers --> showboundaries
    @series l, :sites
    !showbonds && (bonds --> ())
    @series l, :bonds
    showboundaries && @series l, :boundaries
end

const BravaisLatticeValue{Sym, N} = LatticeValue{<:Number, <:OnSites{BravaisLattice{N, <:UnitCell{Sym, N}}}}
raw"""
    rectified_values(lv::LatticeValue)

Returng an array of the values of `lv` on its macrocell.
The $i$-th element of the array corresponds to the $i$-th site of the macrocell.
If the element is `NaN`, it means that the corresponding site is not present in the `lv`'s lattice.

This function might be quite useful in custom plot recipes.
"""
function rectified_values(lv::BravaisLatticeValue{Sym, N} where Sym) where N
    l = stripparams(lattice(lv))
    mins = Vector(l[1].latcoords)
    maxs = Vector(l[1].latcoords)
    for lp in l.pointers
        @. mins = min(mins, lp.latcoords)
        @. maxs = max(maxs, lp.latcoords)
    end
    newvals = fill(NaN, Tuple(@. maxs - mins + 1))
    smins = SVector{N}(mins) .- 1
    for (i, lp) in enumerate(l.pointers)
        newvals[(lp.latcoords - smins)...] = lv.values[i]
    end
    Tuple(mins[i]:maxs[i] for i in 1:N), newvals
end

"""
    plot_fallback(lv::BravaisLatticeValue)

Creates a copy of `lv` lattice value with the underlying `LatticeSym` overwritten to `:GenericBravaisLattice`.
Use it to invoke the default plot recipe for `LatticeValues` when defining a custom one.
"""
function plot_fallback(lv::BravaisLatticeValue)
    l = stripparams(lattice(lv))
    new_l = BravaisLattice(
        UnitCell{:GenericBravaisLattice}(l.unitcell.translations, l.unitcell.basissites),
        l.pointers)
    LatticeValue(new_l, lv.values)
end

@recipe function f(lv::BravaisLatticeValue{:SquareLattice, N}) where N
    N < 3 && (seriestype --> :heatmap)
    uc = unitcell(lattice(lv))
    if plotattributes[:seriestype] === :heatmap
        axes, heatmap_values = rectified_values(lv)
        c_axes = Tuple(axes[i] .+ uc.basissites[i, 1] for i in 1:N)
        aspect_ratio := :equal
        c_axes..., transpose(heatmap_values)
    else
        plot_fallback(lv)
    end
end
