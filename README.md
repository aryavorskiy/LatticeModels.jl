<p align="center"><img src="docs/src/assets/logo.svg" width="250" /></p>

# LatticeModels.jl
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://aryavorskiy.github.io/LatticeModels.jl/dev)
[![codecov](https://codecov.io/gh/aryavorskiy/LatticeModels.jl/branch/master/graph/badge.svg?token=KQN77BL9UV)](https://codecov.io/gh/aryavorskiy/LatticeModels.jl)

This package provides a set of tools to simulate different quantum lattice systems.

WARNING: This package is gradually becoming stable. No new core functional will be introduced. Short list of things to do until final release:
- [ ] Tests
- [ ] Documentation
- [ ] Examples
- [ ] Benchmarks - compare with Kwant and Pybinding

## Installation

Paste the following line into the Julia REPL:
```
]add https://github.com/aryavorskiy/LatticeModels.jl
```
or
```julia
import Pkg; Pkg.add(url="https://github.com/aryavorskiy/LatticeModels.jl")
```

## Similar packages

There are many packages with similar functionality, such as [Quantica.jl](https://github.com/pablosanjose/Quantica.jl), [pybinding](https://docs.pybinding.site/en/stable/index.html) and [Kwant](https://kwant-project.org/). 
However, the scope of these packages is different:

- **Schroedinger equation solvers with time-dependent hamiltonians**. Only Kwant provides similar functionality 
    with its `Tkwant` module, but it lacks performance and flexibility in some cases.
- **Convenient tools for setting boundary conditions and gauge fields**. The only way to do this in 
    Kwant or Pybinding is to manually set the hopping values.
- **A flexible interface for defining new types of lattices and bonds**. Random lattices can be implemented
    on top of `GenericLattice` with ease.
- **Manybody computations**. Kwant and Pybinding are designed mostly for single-particle simulations, while 
    `LatticeModels.jl` can handle manybody systems with particle interaction.

## Example

```julia
using LatticeModels, Plots

# First create a lattice
l = SquareLattice(10, 10)

# Define a tight-binding model hamiltonian with a point flux field through point (5.5, 5.5)
h(B) = tightbinding_hamiltonian(l, field=PointFlux(B, (5.5, 5.5)))

# Find density matrix for filled bands (e. g. with energy < 0)
P_0 = densitymatrix(h(0), mu = 0)

# Perform unitary evolution
τ = 10
a = Animation()
ev = Evolution(t -> h(0.2 * min(t, τ) / τ), P_0)
for state in ev(0:0.1:2τ)
    P, H, t = state
    p = plot(layout=2, size=(800, 400))
    # Find the local density and plot it
    plot!(p[1], localdensity(P), clims=(0, 1), st=:shape, c=:matter)

    # Show currents on the plot
    plot!(p[2], DensityCurrents(H, P), clims=(0, 0.1))

    # Some more tweaks to the plot...
    title!("t = $t")
    frame(a, p)
end

gif(a, "animation.gif")
```

This code creates an animation which displays local density and currents on a plot:

![](animation.gif)