# Examples

## Lowest states in a tight-binding model

The tight-binding model Hamiltonian is defined by this formula:

```math
\hat{H} = \sum_i^\text{sites} \left( c^\dagger_{i + \hat{x}} c_i + c^\dagger_{i + \hat{y}} c_i + h. c. \right)
```

Here we will find its eigenstates and plot their local density on heatmaps.

```@example
using LatticeModels
using Plots
# Generate a 40x40 square lattice
l = SquareLattice(40, 40)
# Define the tight-binding model Hamiltonian
H = tightbinding_hamiltonian(l)

# Calculate eigenvalues and eigenvectors
diag = diagonalize(H)

n = 5
clims = (0, 0.0045)
p = plot(layout = @layout[grid(n, n) a{0.1w}], size=(1000, 850))
for i in 1:n^2
    E_rounded = round(diag.values[i], sigdigits=4)
    plot!(p[i], localdensity(diag[i]), title="\$E_{$i} = $E_rounded\$", st=:shape, 
        clims=clims, c=:inferno, cbar=:none, lw=0, framestyle=:none, xlab="", ylab="")
end

# The following lines are kinda hacky; they draw one colorbar for all heatmaps
plot!(p[n^2+1], framestyle=:none)
scatter!([NaN], zcolor=[NaN], clims=clims, leg=:none, cbar=:right, subplot=n^2+2, 
    background_subplot=:transparent, framestyle=:none, inset=bbox(0.0, 0.05, 0.95, 0.9))
savefig("local_density.png")
nothing # hide
```
![](local_density.png)

## Currents on a ring-shaped sample

In this example we create a ring-shaped sample of a triangular lattice.
Then we adiabatically turn on magnetic field through the hole and see currents emerge.

```@example
using LatticeModels
using Plots

l = TriangularLattice(Circle(10), !Circle(5))
removedangling!(l)
h(B) = tightbinding_hamiltonian(l, field=PointFlux(B))
diag = diagonalize(h(0))

# Find density matrix for filled bands (e. g. energy < 0)
P_0 = densitymatrix(diag, mu = 0)
# Perform unitary evolution
τ = 10
a = Animation()
ev = Evolution(t -> h(0.1 * min(t, τ) / τ), P_0)
for state in ev(0:0.1:2τ)
    P, H, t = state
    # Find the density and plot it
    p = plot(layout=2, size=(1000, 500))
    plot!(p[1], localdensity(P), clims=(0, 1), st=:shape)

    # Show currents on the plot
    plot!(p[2], DensityCurrents(H, P), clims=(0, 0.005), lw=1, arrowheadsize=0.3)

    title!("t = $t")
    frame(a)
end

gif(a, "adiabatic_flux.gif")
```

## Time sequences

In this example we will see how to use `TimeSequence` to store and manipulate time-dependent data.

We will calculate the evolution of a ground state of a tight-binding model after a magnetic field is turned on.
We will store the local density at each time step and use it to plot the local density depending on time, as well as its time derivative and integral over time.

```@example
using LatticeModels, Plots

l = SquareLattice(20, 20)
H = tightbinding_hamiltonian(l)
psi_0 = groundstate(H)
H1 = tightbinding_hamiltonian(l, field=LandauGauge(0.1))
ev = Evolution(H1, psi_0)

densities = TimeSequence{LatticeValue}()
for (psi, _, t) in ev(0:0.1:10)
    densities[t] = localdensity(psi)
end

site_bulk = l[!, x = 10, y = 10]
site_edge = l[!, x = 10, y = 1]
ds_bulk = densities[site_bulk]
ds_edge = densities[site_edge]
plot(ds_bulk, label="ρ(t) (bulk)")
plot!(differentiate(ds_bulk), label="dρ(t)/dt (bulk)")
plot!(ds_edge, label="ρ(t) (edge)")
plot!(integrate(ds_edge), label="∫ρ(t)dt (edge)")
```

## Hofstadter butterfly

The Hofstadter butterfly is a fractal-like structure that appears when the tight-binding model is subjected to a magnetic field. It is a plot of the energy spectrum as a function of the magnetic flux through the unit cell.

To create the Hofstadter butterfly, we will use the Landau gauge for the magnetic field. Note that we have to set periodic boundary conditions, and to make them compatible with the gauge field, they should be tweaked a little:

```math
\psi(x + L_x, y) = \psi(x, y) e^{-2\pi i B y L_x},
\psi(x, y + L_y) = \psi(x, y)
```

Let us plot the Hofstadter butterfiles for square, triangular and honeycomb lattices. The magnetic field field will be changed from zero to one ``\phi_0`` flux quantum per plaquette.

```@example
using LatticeModels, Plots

function get_butterfly(l, lx, ly, plaquette_area)
    xs = Float64[]
    ys = Float64[]
    area = lx * ly
    dflux = 1 / area
    totflux = 1 / plaquette_area
    for B in 0:dflux:totflux
        # magnetic boundary conditions
        f(site) = exp(2pi * im * B * site.y * lx)
        lb = setboundaries(l, [lx, 0] => f, [0, ly] => true)
        H = tightbinding_hamiltonian(lb, field=LandauGauge(B))
        dg = diagonalize(H)
        append!(xs, fill(B, length(dg.values)))
        append!(ys, dg.values)
    end
    return xs, ys
end

p = plot(layout = @layout[a b; _ c{0.5w} _], size=(800, 500), leg=false,
    xlabel="B", ylabel="E")
scatter!(p[1], title="Square lattice",
    get_butterfly(SquareLattice(10, 10), 10, 10, 1), ms=1)
scatter!(p[2], title="Triangluar lattice",
    get_butterfly(TriangularLattice(10, 10), 10, 5 * sqrt(3), sqrt(3) / 4), ms=1)
scatter!(p[3], title="Honeycomb lattice",
    get_butterfly(HoneycombLattice(10, 10), 10, 5 * sqrt(3), sqrt(3) / 2), ms=1)
savefig("hofstadter_butterfly.png")
nothing # hide
```
![](hofstadter_butterfly.png)

## LDOS animation

Local density can be a bit ambiguous for degenerate eigenstates. That's where the LDOS (e. g. the Local Density of States) will be helpful.

The formula for the LDOS is the following:

```math
\text{LDOS}_\alpha(E) = \text{Im} G_{\alpha\alpha}(E - i\delta)
```

where ``G`` is the Green's function and ``\delta`` is the broadening.

Let's create an animation presenting the DOS and LDOS for a square lattice with a hole indside.
We will use the QWZ model hamiltonian, because it has a two-zone band structure, which will make
the results more interesting. See [`qwz`](@ref) for more information about the QWZ model.

```@example
using LatticeModels
using Plots
l = SquareLattice(20, 20)
l_center = l[j1 = 8..13, j2 = 8..13]
setdiff!(l, l_center)   # remove the center
H = qwz(l)

dg = diagonalize(H)
δ = 0.1
Es = -4:0.1:4
Es_d = -4:0.01:4
G = greenfunction(dg)
a = @animate for E in Es
    print("\rE = $E") # hide
    p = plot(layout=2, size=(800, 400))
    plot!(p[1], Es_d, dos(G, broaden=δ), lab="", title="DOS")
    vline!(p[1], [E], lab="")
    plot!(p[2], ldos(G, E, broaden=δ), st=:shape, 
        c=:inferno, clims=(0, NaN), title="LDOS", lw=0)
    plot!(p, plot_title="E = $E, δ = $δ")
end

gif(a, "ldos_animation.gif", fps=10)
```
