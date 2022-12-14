<p align="center"><img src="docs/src/assets/logo.svg" width="250" /></p>

# LatticeModels.jl
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://aryavorskiy.github.io/LatticeModels.jl/)
[![codecov](https://codecov.io/gh/aryavorskiy/LatticeModels.jl/branch/master/graph/badge.svg?token=KQN77BL9UV)](https://codecov.io/gh/aryavorskiy/LatticeModels.jl)

This package provides a set of tools to simulate different quantum lattice systems.

## Installation

```jldoctest
julia> # type ] in REPL to enter pkg mode; then type the following command
pkg> add https://github.com/aryavorskiy/LatticeModels.jl
```
or
```julia
import Pkg; Pkg.add(url="https://github.com/aryavorskiy/LatticeModels.jl")
```

## Sample workflow

```julia
using LatticeModels
using Plots

# First create a lattice
l = SquareLattice(10, 10)

# Define a tight-binding model hamiltonian
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
    plot(site_density(P), clims=(0,1))

    # Show currents on the plot
    plot!(DensityCurrents(H, P), arrows_scale=7)

    # Some more tweaks to the plot...
    print("\rt = $t")
    title!("t = $t")
    frame(a)
end

gif(a, "animation.gif")
```

This code creates an animation which displays local density and currents on a heatmap:

![](animation.gif)