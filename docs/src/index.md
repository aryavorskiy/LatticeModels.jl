# LatticeModels.jl

This package provides a set of tools to simulate different quantum lattice systems.

## Package features
- Bravais lattices with arbitrary geometry and any possible count of internal states on one sites.
- Versatile hamiltonian generation tools.
- Backend-independent computations: linear operators can be of any array type, allowing to use sparse or GPU arrays when needed.
- Smart unitary evolution macro reducing excessive computations where possible.
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) integration.

## Usage examples

### Currents in Hofstadter model on a ring-shaped sample

In this example we delete part of the sites in the middle of a square lattice. 
Then we adiabatically turn on magnetic field through the hole and see currents emerge.

The Hofstadter model hamiltonian is evaluated according to this formula:

$$\hat{H} = \sum_\text{x-links} c^\dagger_i c_j + \sum_\text{y-links} c^\dagger_i c_j + h. c.$$

```@example
using LatticeModels
using LinearAlgebra, Plots

l = SquareLattice(10, 10) do site, (x, y)
    abs(x) > 1 || abs(y) > 1
end

# Define a Hofstadter model hamiltonian
h(B) = @hamiltonian begin   
    lattice := l
    # Add hoppings along axis x and y
    @hop axis = 1
    @hop axis = 2
    # Add magnetic field through (0, 0) point
    field := FluxField(B, (0, 0))
end

# Calculate eigenvalues and eigenvectors
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
    density = diag_aggregate(m -> real(tr(m)), P)
    heatmap(density, clims=(0,1))

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
\sum_\text{x-links} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_\text{y-links} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
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
    p = plot(layout=2, size=(900, 500))

    # Local Chern marker heatmap
    lcm_operator = 4pi * im * P * X * P * Y * P
    chern_marker = diag_aggregate(tr, lcm_operator) .|> real
    heatmap!(p[1], chern_marker, clims=(-2, 2))

    # Select sites on y=0 line (use ≈ to avoid rounding errors)
    chern_marker_on_sw = chern_marker[@. y ≈ 0]
    # Mark selected sites on the heatmap
    plot!(p[1], chern_marker_on_sw.lattice, high_contrast=true)
    # Add a line plot
    plot!(p[2], chern_marker_on_sw, project_axis=:x, ylims=(-3, 3), lab=:none)

    plot!(plot_title="t = $t")
    frame(a)
end

gif(a, "animation.gif")
```