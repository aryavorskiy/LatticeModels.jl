# Lattices

This chapter describes the basic functionality of `LatticeModels.jl` and how to use it to create and manipulate lattices.

## Basics

Creating a lattice is simple. If you want to work with Bravais lattices, there are plenty of predefined types in the package. For example, to create a square lattice, you can use the `SquareLattice` type:

```@example 1
using LatticeModels, Plots
l = SquareLattice(10, 10)
plot(l)
```

This will create a 10x10 square lattice. Note that you can create a plot by simply using the `plot` function from the `Plots` package.

`SquareLattice(10, 10)` notation means that we take a square lattice unit cell and translate it 10 times in the x and y directions. This syntax can be extended a little bit:

```@example 1
l = HoneycombLattice(-2:2, -2:2)
plot(l)
plot!(l[j1 = 0, j2 = 0], c=:red, ms=6)
```

What happened here? We created a honeycomb lattice by translating the unit cell from -2 to 2 in both directions. Then we plotted the lattice and highlighted the unit cell at the origin to make it more visible. `l[j1 = 0, j2 = 0]` allowed us selecting part of the lattice by specifying the indices of the unit cell (`j1` and `j2` are the indices of the unit cell in the first and second directions, respectively).

If we need to create a lattice with a less trivial shape, we can use any function we need:

```@example 1
l = TriangularLattice(-10:10, -10:10) do site
    return 4 < sqrt(site.x^2 + site.y^2) < 8     # Create a ring-shaped lattice
end
plot(l)
```

Here the lattice constructor first translated the unit cell from -10 to 10 in both directions, and then applied the function to each site to create a ring-shaped lattice. This is similar to the `filter` function in Julia — in fact, you can use `filter` or `filter!` on an existing lattice to create a new one as well.

There are other things you can control when creating a lattice, such as the lattice offset & rotation.

```@example 1
l1 = SquareLattice(-2:2, -2:2)
l2 = SquareLattice(-2:2, -2:2, offset=:centeralign)
l3 = SquareLattice(-2:2, -2:2, rotate=pi/3, offset=[6, -1.5])
plot(l1, lab="No offset", shape=:circle)
plot!(l2, lab="Center the unit cell", shape=:star)
plot!(l3, lab="Shifted and rotated by π/3", shape=:square)
```

To find out more about offset and rotation, see [`UnitCell`](@ref) — the keywords are described there.

!!! note
    If you use both offset/rotation and a function to create a lattice, the function will be applied to the sites **after** the offset/rotation is applied. Use the `postoffset` and `postrotate` keywords to control the position and orientation of the lattice after the function is applied.
    
The lattices implement the `AbstractSet` interface, so you can use all the set operations on them — `union`, `intersect`, `setdiff` etc.

```@example 1
l1 = SquareLattice(-2:0, -2:0)
l2 = SquareLattice(0:2, 0:2)
l3 = SquareLattice(-3:3, -3:3)
l = setdiff(l3, union(l1, l2))
plot(l)
plot!(union(l1, l2), showbonds=false, alpha=0.3)
```

## Sites

Let's find out what sites actually are. A site is generally a point in the lattice. It is defined by its position in space and maybe some additional properties. In case of a Bravais lattice these additional properties are unit cell indices and the index _in_ the unit cell.

A lattice is generally a set-like structure that allows indexing. Let's take a closer look in the REPL:

```@repl 2
using LatticeModels, Plots
l = HoneycombLattice(-2:2, -2:2)
site = l[1]     # Get the first site
site.x          # Get the x-coordinate of the site
site.j2         # Get the second index of the unit cell
site.index      # Get the index of the site in the unit cell
x, y = site     # Destructure the site
```

As we see, we can access the properties of the site simply as fields of the site object. We can also destructure the site to get its coordinates.

!!! compat "Julia 1.8"
    Accessing the properties of the site as fields like `site.x` requires Julia 1.8 and later. This limitation is imposed with purpose, since this seriously affects runtime performance in earlier versions. You will still be able to destructure the site to get its coordinates, or use the following fields:
    - `site.coords` — the position of the site
    - `site.latcoords` — the unit cell indices
    - `site.index` — the index of the site in the unit cell

Properties like `x`, `j1`, `index` etc. are part of a general `SiteProperty` interface. You can use them to create 'slices' of lattices:

```@example 2
slice = l[j1 = 0..2, j2 = -2..0, index=1]  # Get a slice of the lattice
plot(l)
plot!(slice, c=:red, ms=6)
```

Here `0..2` and `-2..0` are intervals defining the ranges of the unit cell indices. You can use any collection instead of then if you need.

Finding sites by their properties can be done with the same notation:

```@repl 2
l[x = 1.5, y = √3/2]    # Find the site with x = 1.5 and y = √3/2
l[x = 1.2, y = 3]       # No such site, throws an error
```

This is notation is convenient yet type-unstable, since it returns a `Site` object if there is one site satisfying the condition — otherwise a lattice is returned. To make sure that the result is indeed a site, add `!` to the beginning of the condition:

```@repl 2
l[!, x = 1.5, y = √3/2]    # Find the site with x = 1.5 and y = √3/2
l[!, x = 1.5]              # More than one site, throws an error
```

## Custom `UnitCell`

You can also create a lattice from a custom unit cell:

```@example 3
using LatticeModels, Plots
# This will be our custom honeycomb lattice unit cell
# First argument - vectors of the unit cell
# Second argument - radius-vectors for the sites in the unit cell
uc = UnitCell([[1/2, sqrt(3)/2] [-1/2, sqrt(3)/2]], [[0, sqrt(3)/6] [0, -sqrt(3)/6]])
plot(uc)    # Plot the unit cell
```

Note that both arguments are actually matrices — the first one is a matrix of the unit cell vectors, and the second one is a matrix of the site positions in the unit cell. However, here we used concatenation to create the matrices for the sake of readability: remember that `[[a, b] [c, d]]` is equivalent to `[a c; b d]`.

To create a lattice, we can use the [`span_unitcells`](@ref) function:

```@example 3
l = span_unitcells(uc, -5:5, -5:5) do site
    x, y = site
    return abs(y) < 5 && 
        abs(y * 1 / 2 + x * sqrt(3) / 2) < 5 && 
        abs(y * 1 / 2 - x * sqrt(3) / 2) < 5
end   # Create a hex shape
plot(l)
```

In fact, the constructors we discussed earlier are just a shorthand for `span_unitcells` with a predefined unit cell. You can use `span_unitcells` to create a lattice from any unit cell you want.

## Shapes

The shapes framework is a powerful tool for creating lattices of arbitrary geometry:

```@example 4
using LatticeModels, Plots
l = SquareLattice{2}(Hexagon(10, [-10, 0]), Circle(10, [10, 0]))
plot(l)
```

Here we created a square lattice in shape of a hexagon and a circle. The first argument of the shape is its radius (for the hexagon it is the distance from its center to the vortices), and the second argument is the center of the shape. Other possible shapes include `Rectangle`, `Polygon`, `SiteAt` and `Path`.

```@example 4
complex_l = SquareLattice{2}(   # Here you have to specify the dimension of the lattice
    Circle(10), Circle(10, [20, 0]), Circle(10, [10, 10√3]),
    !Circle(5), !Circle(5, [20, 0]), !Circle(5, [10, 10√3]),
    Rectangle(-5..5, -14..(-12)), Rectangle(15..25, -14..(-12)),
    Path([-12, 32], [32, 32]), SiteAt([0, 0]), SiteAt([20, 0])
)
plot(complex_l)
```

Note that adding `!` before the shape inverts it. This is useful when you need to create a lattice with a hole in it.

Sometimes the shape can become ill-formed — this happens when the unit cell has non-trivial geometry. In this case you may need to remove the dangling sites using the [`removedangling!`](@ref) function. For example, they can arise when creating a path on a honeycomb lattice:

```@example 4
l = HoneycombLattice(Circle(3, [0, 0]), Circle(3, [-2, 10]), Path([0, 0], [-2, 10]))
p = plot(size=(800, 350), layout=(1, 2))
plot!(p[1], l, title="With dangling sites")
removedangling!(l)
plot!(p[2], l, title="Without dangling sites")
```

Let's discuss what is happening under the hood. The `HoneycombLattice` constructor calls the [`fillshapes`](@ref) function, which estimates the unitcells one has to span for each shape, and adds the sites that are in the shape to the lattice. 

You can also use `addshapes!` to add shapes to an existing lattice and `deleteshapes!` to remove them. These functions, however, do not support the `!` notation for inverting shapes.

One last, but not least, thing to mention is that this framework allows approximate scaling. If you need a lattice with distinct shape and, say, roughly 1000 sites, you can use the `sites` keyword to specify the number of sites you need:

```@example 4
l = TriangularLattice(sites=1000, Circle(1, [-1, 0]), Circle(1, [1, 0]), Circle(1, [0, √3]))
plot(l, title = "$(length(l)) ≈ 1000 sites")
```

You can also use the [`shaperadius`](@ref) function to estimate the radius of the shape that will give you the desired number of sites:[^1]

```@example 4
circ = SquareLattice{2}(Circle(), sites=150)
r = shaperadius(circ, Circle())
plot(circ, lab = "$(length(circ)) ≈ 150 sites")
plot!(Circle(r), c=:grey, ls=:dash, lab = "r ≈ $r")
```

[^1]: Actually, the `shaperadius` function returns the scaling factor for the shape set. However, `Circle()` by default creates a circle with radius 1, so the scaling factor is equal to the radius of the circle.

!!! note
    Radius estimation is not always precise and works under following assumptions:
    - The shapes are large enough to contain the unit cell and do not intersect with each other.
    - The inverted shapes are all contained in the non-inverted ones, and also do not intersect with each other.

## Bonds and hoppings

In the scope of this package, bonds are not considered part of the lattice, but rather a separate structure that connects sites. You can assign nearest-neighbour hoppings to the lattice (of course, since you can see them on the plot), but in general the lattice and the bonds are separate entities. The reasoning behind this is that in some cases you may need to use different sets of bonds for different models.

To find out more about bonds, adjacency and boundary conditions, see the next chapter: [Adjacency and boundary conditions](@ref).
