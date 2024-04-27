# Currents

This chapter covers the calculation of currents in the system.

## Basics

An `AbstractCurrents` object is a collection of currents of some nature (e.g. charge, spin, etc.) from one site to another. 
In most cases it is a lazy object that doesn't store the currents themselves, but rather the parameters to evaluate them on the go.

A good example of this is [`DensityCurrents`](@ref), which calculates the charge currents between sites.

```@example 1
using LatticeModels, Plots
l = SquareLattice(10, 10)
H = tightbinding_hamiltonian(l)
P = densitymatrix(H, statistics=FermiDirac, mu=0)
H1 = tightbinding_hamiltonian(l, field=PointFlux(0.1, (5.5, 5.5)))
C = DensityCurrents(H1, P)
plot(C)
```

Here we find the density matrix `P` of the system, add a magnetic field to the Hamiltonian to induce currents (`H1`), and calculate the charge currents between sites (`C`). The resulting plot shows the currents between sites.

Currents can be indexed with sites to get the currents between them:

```@repl 1
site1 = l[!, x=1, y=1]; site2 = l[!, x=1, y=2]; site3 = l[!, x=1, y=3];
C[site1, site2]
C[site1, site3]         # Current between site1 and site3 is zero
```

You can evaluate all currents at once by converting the currents object to [`Currents`](@ref). This can be useful when you want to perform non-trivial operations on the currents:

```@repl 1
C1 = Currents(C)
C1[site1, site2] == C[site1, site2]
```

Also note that the `Currents` object is writable:

```@example 1
site4 = l[!, x=2, y=2]
C1[site1, site2] = 0.05
C1[l[!, x=10, y=10], l[!, x=9, y=9]] = 0.05
plot(C1)
```

## Common Operations

One function you may find useful is [`currentsfromto`](@ref), which calculates currents between two sites or other domains.
As an example, let's find how the currents through a circular sample depends on time:

```@example 2
using LatticeModels, Plots
l = SquareLattice{2}(Circle(10), !Circle(5))
h(t) = tightbinding_hamiltonian(l, field=PointFlux(0.1 * t, (0, 0)))
P = densitymatrix(h(0), statistics=FermiDirac, mu=0)
τ = 10
ev = Evolution(t -> h(min(t, τ)), P)
dom1 = l[x = -Inf .. 0, y = 0]
dom2 = l[x = -Inf .. 0, y = 1]
currs = TimeSequence(ev, 0:0.1:2τ) do (Pt, Ht, t)
    Currents(DensityCurrents(Ht, Pt))    # Convert to Currents, since P is reused
end
currs_from_to = map(curr -> currentsfromto(curr, dom1, dom2), currs)
plot(currs_from_to, title="Currents through the circle")
```

Here we perform a time evolution of the system with a magnetic field that is slowly increasing from 0 to `0.1 * τ`, after which it remains constant. See the [Evolution](@ref evolution_chapter) chapter for more details about time evolution.

The value we evaluated is the currents between two domains `dom1` and `dom2` at each time step. Here are 
the domains:

```@example 2
C = currs[τ]
plot(C)
plot!(dom1, label="Domain 1")
plot!(dom2, label="Domain 2")
```

Another useful function is [`currentsfrom`](@ref), which calculates currents from a given domain (site or part of the lattice) to the rest of the lattice and returns a [`LatticeValue`](@ref) object.

Let us continue with the previous example. This will calculate the currents from `site1` in the bottom-left corner to the rest of the lattice:

```@example 2
site = l[!, x = -7, y = 2]
heatmap(currentsfrom(C, site), c=:balance, clims=(-0.02, 0.02))
plot!(site, label="", title="Currents from site $(site.coords)")
```

Also note that currents are iterable objects, yielding a tuple of `source => target` pair and the current value. 
The order of the sites is such that the current is directed from the source to the target and is always positive.
This can be useful to analyze, for example, the direction of the currents:

```@example 2
cs = zeros(4)
ns = zeros(Int, 4)
for ((src, tgt), curr) in C
    @assert curr >= 0
    curr < 1e-6 && continue
    x, y = normalize(tgt.coords - src.coords)
    angle = atan(y, x)
    idx = mod(round(Int, angle / (π / 2)), 4) + 1
    cs[idx] += curr
    ns[idx] += 1
end
bar(["0", "π/2", "π", "3π/2"], cs ./ ns, label="Average current by direction")
```

## Visualization

The currents can be visualized using the `plot` function. The values of the currents will be shown as arrows between the sites. The color of the arrows represents the magnitude of the currents. Also note that the opacity of the arrows is also proportional to the magnitude of the currents, which can be disabled by setting `arrowtransparency=false`. You can customize the appearance of the arrows using the `arrowheadsize` and `arrowheadwidth` keywords.

Here we will continue with the `C` currents obtained in the previous example:

```@example 2
p = plot(layout = @layout[a b; _ c{0.5w} _], size=(1000, 1000))
plot!(p[1], C, title="Default")
plot!(p[2], C, arrowtransparency=false, title="Opaque arrows")  
plot!(p[3], C, arrowheadsize=0.3, arrowheadwidth=1, arrowtransparency=false,
    title="Custom arrowhead")
```

You can plot an `AbstractCurrents` and a `LatticeValue` object in one axes. However, due to Plots.jl limitations, both series must share the same colorbar. You can bypass this by plotting them on separate axes, or by plotting the currents with a single color. In this case you need to set `line_z=nothing` to disable the colorbar for the arrows. Here is an example of plotting the currents and the lattice together:

```@example 2
plot(localdensity(P))
plot!(C, c=:lawngreen, arrowheadsize=0.3, arrowheadwidth=1, line_z=nothing)
```
