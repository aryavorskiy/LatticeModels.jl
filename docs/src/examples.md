# Examples

## Local density for lowest states in a tight-binding model

The tight-binding model hamiltonian is defined by this formula:

$$\hat{H} = \sum_i^\text{sites} \left( c^\dagger_{i + \hat{x}} c_i + c^\dagger_{i + \hat{y}} c_i + h. c. \right)$$

Here we will find its eigenstates and plot their local density on heatmaps.

```@example
using LatticeModels
using Plots
# Generate a 40x40 square lattice
l = SquareLattice(40, 40)
# Define the tight-binding model hamiltonian
H = tightbinding_hamiltonian(l)

# Calculate eigenvalues and eigenvectors
diag = diagonalize(H)

n = 5
clims = (0, 0.0045)
p = plot(layout = @layout[ grid(n, n) a{0.1w}], size=(1000, 850))
for i in 1:n^2
    E_rounded = round(diag.values[i], sigdigits=4)
    plot!(p[i], localdensity(diag[i]), title="\$E_{$i} = $E_rounded\$", clims=clims, cbar=:none)
end

# The following 2 lines are kinda hacky; they draw one colorbar for all heatmaps
plot!(p[n^2+1], framestyle=:none)
scatter!([NaN], zcolor=[NaN], clims=clims, leg=:none, cbar=:right, background_subplot=:transparent, 
    framestyle=:none, inset=bbox(0.0, 0.05, 0.95, 0.9), subplot=n^2+2, c=:matter)
```

## Currents in a tight-binding model on a ring-shaped sample

In this example we create a ring-shaped sample of a triangular lattice.
Then we adiabatically turn on magnetic field through the hole and see currents emerge.

The tight-binding hamiltonian is the same as in the example above.

```@example
using LatticeModels
using Plots

l = TriangularLattice(Circle(10), !Circle(5))
removedangling!(l)
h(B) = tightbinding_hamiltonian(l, field=FluxField(B))
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
    p = plot(layout=2, size=(800, 400))
    plot!(p[1], localdensity(P), clims=(0, 0.1), st=:shape)

    # Show currents on the plot
    plot!(p[2], DensityCurrents(H, P))

    title!("t = $t")
    frame(a)
end

gif(a, "animation.gif")
```

## Local Chern marker with hamiltonian quench

The Chern insulator hamiltonian is described by this formula:

$$\hat{H} = 
\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i + 
\sum_i^\text{sites} \left( 
c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i + 
c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i + 
h. c. \right)$$

In this experiment we create a filled state density matrix for a system with $m_i = 1$. 
After that we change the $m_i$ in the center of the lattice to $-1$ and start the evolution.

```@example
using LatticeModels
using Plots

l = SquareLattice(11, 11)
x, y = coordvalues(l)

# Initial hamiltonian: m=1 everywhere
H1 = qwz(l, 1)

# Quenched hamiltonian: m=-1 in the central 3x3 square
M = ones(l)
M[x = 4..8, y = 4..8] .= -1
H2 = qwz(M)
X, Y = coordoperators(l, 2)

sp = diagonalize(H1)
P_0 = densitymatrix(sp, mu = 0)

τ = 10
a = Animation()
ev = Evolution(H2, P_0)
for state in ev(0:0.1:2τ)
    P, H, t = state
    p = plot(layout=2, size=(900, 400))

    # Local Chern marker heatmap
    lcm_operator = 4pi * im * P * X * P * Y * P
    chern_marker = localdensity(lcm_operator)
    plot!(p[1], chern_marker, clims=(-2, 2), st=:shape)

    # Select sites on y=6 line
    chern_marker_on_sw = chern_marker[y = 6]
    # Mark selected sites on the plot
    plot!(p[1], lattice(chern_marker_on_sw), high_contrast=true)
    # Add a line plot
    plot!(p[2], project(chern_marker_on_sw, :x), ylims=(-3, 3), lab=:none)

    plot!(plot_title="t = $t")
    frame(a)
end

gif(a, "animation.gif")
```

## LDOS animation

Local density can be a bit ambiguous for degenerate eigenstates. That's where the LDOS (Refer to [`ldos`](@ref) documentation) will be helpful.

Let's take the same hamiltonian from the previous example and create a LDOS animation.

```@example
using LatticeModels
using Plots
l = SquareLattice(40, 40)
H = qwz(l, 1)

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
    plot!(p[2], ldos(G, E, broaden=δ), clims=(0, NaN), title="LDOS")
    plot!(p, plot_title="E = $E, δ = $δ")
end

gif(a, "animation.gif", fps=10)
```