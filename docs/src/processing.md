This chapter is devoted to various ways of working with values defined on lattice sites. The [`LatticeValue`](@ref) is 
a struct designed specially for this.

## Generate lattice-dependent values

You can create new `LatticeValue`s mostly the same way as with ordinary `Vector`s:

```@repl env
using LatticeModels, Plots
l = SquareLattice(10, 10)
lv = zeros(l)                       # 0.0 on all sites
lv1 = ones(l)                       # 1.0 on all sites
lv2 = rand(l)                       # uniformly distributed random numbers
lv3 = randn(l)                      # normally distributed random numbers
nothing # hide
```

To generate a tuple of `LatticeValue`s for site coordinates, you can use the [`coord_values`](@ref) function. Also use `param_value` to generate `LatticeValue`s representing various site parameters similar to [`param_operator`](@ref):

```@repl env
x, y = coord_values(l)
x1 = param_value(l, :x)     # Same as `x`
j1 = param_value(l, p"j1")  # Unit cell index #1
scatter(j1)                 # Show values on a scatter plot
```

Another way to generate a `LatticeValue` is using broadcasting or `do`-syndax. These are the most common ways to generate arbitrary `LatticeValue`s that can be used for different models or wavefunctions:

```@repl env
xp2y = @. x + 2y
xp2y = l .|> (site -> site.x + 2 * site.y)
xo2y = LatticeValue((site -> site.x + 2 * site.y), l)
xp2y = LatticeValue(l) do (x, y)
    return x + 2y
end
plot(xp2y)
```

All these examples in this piece of code produce the same `LatticeValue`. In the last example `(x, y)` is the 
parameter of the `do`-lambda, which unpacks the `LatticeSite` in-place.

Note that the plot displaying the `LatticeValue` values is now a heatmap. This is because the default plot recipes 
that is triggered when `plot` is called is overridden for square lattices. However, any lattice type implements a 
scatter plot recipe, and for most lattices this is the default. Consider the honeycomb lattice: 

```@repl env
hl = HoneycombLattice(5, 5)
x, y = coord_values(hl)
plot(@. x + 2y)             # heatmap(@. x + 2y) will produce an error
```

### Indexing

Indexing a `LatticeValue` works the same way as indexing a `Lattice`; it produces a new `LatticeValue` with the same 
values, but defined on a sublattice. If the resulting sublattice consists of only one site, the value on this site is 
returned. Storing new values using the same indexing syntax is also available.

```@example env
l = SquareLattice(10, 10)
x, y = coord_values(l)
lv = zeros(l)
lv[x = 1..3, y = 2..9] = x      # Assign another LatticeValue
lv[x = 7..9, y = 2..9] .= 3     # or a number
nothing # hide
```

Now let us plot a sectional view of this `LatticeValue` along the `y = 5` line. We will need the `project` function to 
convert a `LatticeValue` defined on a "stripe" of sites into a series that can be shown on a line plot:

!!! info
    `project` accepts any site parameter as second argument and returns a tuple of two lists, one containing values of 
    the parameter on the sites of the `LatticeValue`, and the other values of the latter.

```@example env
p = plot(layout=(2, 1), size=(500, 900))
plot!(p[1], lv)
plot!(p[1], l[x = 5], :high_contrast)   # Mark the sites that will be included into the sectional view
plot!(p[2], project(l[x = 5], :y))      # The resulting narrowed `LatticeValue` is projected onto the y axis
```

Another way is to index a `Lattice`/`LatticeValue` using a boolean-typed `LatticeValue`. The behavior is pretty much the same:

```@example env
# Let's select the sites in a sector-shaped area
plot(lv[@. √(x^2 + y^2) ≥ 10])
```

The approach from above provides a flexible way to edit `LatticeValue`s:

```@example env
lv2 = ones(l)
lv2[@. x < y] = x .* y          
lv2[@. x ≥ y && x ≥ -y] .= 20 
plot(lv2)
```

## Generate a wavefunction

!!! info "Docs under construction"
    This manual page is not written yet. I'm working on that. You can watch [a video on YouTube](https://www.youtube.com/watch?v=dQw4w9WgXcQ) meanwhile.
## Obtain lattice density

!!! info "Docs under construction"
    This manual page is not written yet. I'm working on that. You can watch [a video on YouTube](https://www.youtube.com/watch?v=dQw4w9WgXcQ) meanwhile.