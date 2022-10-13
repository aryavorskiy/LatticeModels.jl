# LatticeModels.jl

This package provides a set of tools to simulate different quantum lattice systems.

## Package features
- Bravais lattices with arbitrary geometry and any possible count of internal states on one sites.
- Versatile hamiltonian generation tools.
- Backend-independent computations: linear operators can be of any array type, allowing to use sparse or GPU arrays when needed.
- Smart unitary evolution macro reducing excessive computations where possible.
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) integration.

## Example usage

```@example
using LatticeModels
using LinearAlgebra, Plots

l = SquareLattice(10, 10)

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
    plot(density, clims=(0,1))

    # Show currents on the plot
    plot!(DensityCurrents(H, P), arrows_scale=7)

    title!("t = $t")
    frame(a)
end

gif(a, "animation.gif")
```