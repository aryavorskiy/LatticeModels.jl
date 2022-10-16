## LatticeValue basics

A LatticeValue is a struct that maps sites of a certain `Lattice` to values of some type. 
One can be generated using a `do`-syntax similar to one in [Sublattices](@ref):

```@setup env
using LatticeModels, Plots
```

```@repl env
l = SquareLattice(5, 5)
lv = LatticeValue(l) do site, (x, y); x + y + 1; end
```

To generate `LatticeValue`s for site coordinates, you can use the [`coord_values`](@ref) function.
Note that `LatticeValue`s support broadcasting, which means you can create coordinate-dependent lattice values in-place:

```@repl env
x, y = coord_values(l)
lv == x .+ y .+ 1
```

Lattice values implement a scatter plot recipe, which colors the plot markers according to the value:
```@example env
plot(lv, markersize=20)
```

!!! note
    A wave function cannot and must not be stored as a `LatticeValue` - use `LatticeVector` instead. 
    The reason is that a `LatticeValues` does not take one-site phase spaces into account.

Depending on the lattice type, additional plot recipes can be available. For example, a lattice value on a square lattice can be plotted as a heatmap:

```@example env
heatmap(lv)
```

## Indexing

It is often required to select some sites by certain condition. 
This can be done using a `LatticeValue{Bool}` and broadcasting (like with [Sublattices](@ref)).

```@example env
heatmap(lv[@. √(x^2 + y^2) > 1.2])
```

Note that a `LatticeValue` can be projected to some coordinate axis to create line plots.

```@example env
lv_on_line = lv[@. x ≈ 0]   # Use approximate comparison to avoid rounding errors
p = plot(layout=(2, 1))

heatmap!(p[1], lv)
plot!(p[1], lattice(lv_on_line), high_contrast=true)
plot!(p[2], lv_on_line, project_axis=:y)
```

Note that we can show the sites we selected by plotting the lattice of the selected values with `high_contrast=true`.
This options make the plot markers black-and-white, which prevents them from blending in with the heatmap in the background.