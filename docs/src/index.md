# LatticeModels.jl

This package provides a set of tools to simulate different quantum lattice systems.

## Installation

```jldoctest
pkg> add https://github.com/aryavorskiy/LatticeModels.jl
```
or
```julia
import Pkg; Pkg.install(url="https://github.com/aryavorskiy/LatticeModels.jl")
```

## Package features
- Bravais lattices with arbitrary geometry and any possible count of internal states on one sites.
- Versatile hamiltonian generation tools.
- Backend-independent computations: linear operators can be of any array type, allowing to use sparse or GPU arrays when needed.
- Smart unitary evolution macro reducing excessive computations where possible.
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) integration.

## Usage examples

### Local density for lowest states in a tight-binding model

The tight-binding model hamiltonian is defined by this formula:

$$\hat{H} = \sum_\text{x-bonds} c^\dagger_i c_j + \sum_\text{y-bonds} c^\dagger_i c_j + h. c.$$

Here we will find its eigenstates and plot their local density on heatmaps.

```@example
using LatticeModels
using LinearAlgebra, Plots
# Generate a 40x40 square lattice
l = SquareLattice(40, 40)
# Define the tight-binding model hamiltonian
H = @hamiltonian begin 
    lattice := l
    # Add hoppings along axis x and y
    @hop axis = 1
    @hop axis = 2
end

# Calculate eigenvalues and eigenvectors
sp = spectrum(H)

n = 5
clims = (0, 0.0045)
p = plot(layout = @layout[ grid(n, n) a{0.1w}], size=(1000, 850))
for i in 1:n^2
    E_rounded = round(eigvals(sp)[i], sigdigits=4)
    plot!(p[i], site_density(sp[i]), title="\$E_{$i} = $E_rounded\$", clims=clims, cbar=:none)
end

# The following 2 lines are kinda hacky; they draw one colorbar for all heatmaps
plot!(p[n^2+1], framestyle=:none)
scatter!([NaN], zcolor=[NaN], clims=clims, leg=:none, cbar=:right, background=:transparent, 
    framestyle=:none, inset=bbox(0.0, 0.05, 0.95, 0.9), subplot=n^2+2)
```

### Currents in a tight-binding model on a ring-shaped sample

In this example we delete part of the sites in the middle of a square lattice. 
Then we adiabatically turn on magnetic field through the hole and see currents emerge.

The tight-binding hamiltonian is the same as in the example above.

```@example
using LatticeModels
using LinearAlgebra, Plots

l = SquareLattice(10, 10) do site, (x, y)
    abs(x) > 1 || abs(y) > 1
end
h(B) = @hamiltonian begin   
    lattice := l
    @hop axis = 1
    @hop axis = 2
    # Add magnetic field through (0, 0) point
    field := FluxField(B, (0, 0))
end
sp = spectrum(h(0))

# Find density matrix for filled bands (e. g. energy < 0)
P_0 = filled_projector(sp)
# Perform unitary evolution
τ = 10
a = Animation()
@evolution {
    H := h(0.1 * min(t, τ) / τ)
    P_0 --> H --> P
} for t in 0:0.1:2τ
    # Find the partial trace and plot it
    plot(site_density(P), clims=(0,1))

    # Show currents on the plot
    plot!(DensityCurrents(H, P), arrows_scale=7, arrows_rtol=0.1)

    title!("t = $t")
    frame(a)
end

gif(a, "animation.gif")
```

### Local Chern marker with hamiltonian quench

The Chern insulator hamiltonian is described by this formula:

$$\hat{H} = 
\sum_i m_i c^\dagger_i \sigma_z c_i + \left(
\sum_\text{x-bonds} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_\text{y-bonds} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c. \right)$$

In this experiment we create a filled state density matrix for a system with $m_i = 1$. 
After that we change the $m_i$ in the center of the lattice to $-1$ and start the evolution.

```@example
using LatticeModels
using LinearAlgebra, Plots

l = SquareLattice(11, 11)
x, y = coord_values(l)

# The Pauli matrices
σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]

# Initial hamiltonian: m=1 everywhere
H1 = @hamiltonian begin   
    lattice := l
    @diag σ[3]
    @hop (σ[3] - im * σ[1]) / 2 axis = 1
    @hop (σ[3] - im * σ[2]) / 2 axis = 2
end

# Quenched hamiltonian: m=-1 in the central 3x3 square
M = @. (abs(x) < 1.5 && abs(y) < 1.5) * -2 + 1
H2 = @hamiltonian begin
    lattice := l
    @diag M ⊗ σ[3]
    @hop (σ[3] - im * σ[1]) / 2 axis = 1
    @hop (σ[3] - im * σ[2]) / 2 axis = 2
end
X, Y = coord_operators(l, 2)

sp = spectrum(H1)
P_0 = filled_projector(sp)

τ = 10
a = Animation()
@evolution {
    H := H2
    P_0 --> H --> P
} for t in 0:0.1:2τ
    p = plot(layout=2, size=(900, 400))

    # Local Chern marker heatmap
    lcm_operator = 4pi * im * P * X * P * Y * P
    chern_marker = ptrace(lcm_operator) .|> real
    plot!(p[1], chern_marker, clims=(-2, 2))

    # Select sites on y=0 line (use ≈ to avoid rounding errors)
    chern_marker_on_sw = chern_marker[@. y ≈ 0]
    # Mark selected sites on the heatmap
    plot!(p[1], lattice(chern_marker_on_sw), high_contrast=true)
    # Add a line plot
    plot!(p[2], project(chern_marker_on_sw, :x), ylims=(-3, 3), lab=:none)

    plot!(plot_title="t = $t")
    frame(a)
end

gif(a, "animation.gif")
```

### LDOS animation

Local density can be a bit ambiguous for degenerate eigenstates. That's where the LDOS (Refer to [`ldos`](@ref) documentation) will be helpful.

Let's take the same hamiltonian from the previous example and create a LDOS animation.

```@example
using LatticeModels
using LinearAlgebra, Plots
l = SquareLattice(40, 40)
σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
H = @hamiltonian begin   
    lattice := l
    @diag σ[3]
    @hop (σ[3] - im * σ[1]) / 2 axis = 1
    @hop (σ[3] - im * σ[2]) / 2 axis = 2
end

sp = spectrum(H)
δ = 0.1
Es = -4:0.1:4
Es_d = -4:0.01:4
ldosf = ldos(sp, δ)
a = @animate for E in Es
    print("\rE = $E") # hide
    p = plot(layout=2, size=(800, 400))
    plot!(p[1], Es_d, dos(sp, δ), lab="", title="DOS")
    vline!(p[1], [E], lab="")
    plot!(p[2], ldosf(E), clims=(0, NaN), title="LDOS")
    plot!(p, plot_title="E = $E, δ = $δ")
end

gif(a, "animation.gif", fps=10)
```
