using RecipesBase

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
