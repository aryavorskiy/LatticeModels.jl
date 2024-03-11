# LatticeModels.jl

This package provides a set of tools to simulate different quantum lattice systems.

## Installation

Paste the following line into the Julia REPL:
```
]add https://github.com/aryavorskiy/LatticeModels.jl
```
or
```julia
import Pkg; Pkg.add(url="https://github.com/aryavorskiy/LatticeModels.jl")
```

!!! compat
    This package has tested compatibility with Julia v1.6. Use with caution on lower versions.

## Package features
- Bravais lattices with arbitrary geometry or boundary conditions.
- Powerful operator generation tools.
- Backend-independent computations: linear operators can be of any array type, allowing to use sparse or GPU arrays when needed.
- Manybody computations.
- Smart unitary evolution reducing excessive computations where possible.
- Supports visualization with [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- Compatible with [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Similar packages

Packages such as [Quantica.jl](https://github.com/pablosanjose/Quantica.jl), [pybinding](https://docs.pybinding.site/en/stable/index.html) and [Kwant](https://kwant-project.org/) provide similar functionality. However, they lack some features that are present in `LatticeModels.jl`:

- **Schroedinger equation solvers with time-dependent hamiltonians**. Only Kwant provides similar 
    functionality with its `Tkwant` module, but it lacks the performance and flexibility of the 
    [`Evolution`](@ref) struct.
- **Convenient tools for setting boundary conditions and gauge fields**. The only way to do this in 
    Kwant or Pybinding is to manually set the hopping values.
- **A flexible interface for defining new types of lattices and bonds**. Random lattices can be implemented
    on top of [`GenericLattice`](@ref) with ease.
- **Manybody computations**. Kwant and Pybinding are designed mostly for single-particle simulations, while 
    `LatticeModels.jl` can handle manybody systems with particle interaction.

Still, these packages have their own gigantic advantage - they support leads. This is a feature that is not present in `LatticeModels.jl`, and it is not planned to be implemented in the near future.

## Usage example

This simple code plots local density for lowest energy states of a square tight-binding lattice.

```@example
using LatticeModels, Plots

l = SquareLattice(40, 40)
H = tightbinding_hamiltonian(l)
diag = diagonalize(H)

n = 5
clims = (0, 0.0045)
p = plot(layout = @layout[ grid(n, n) a{0.1w}], size=(1000, 850))
for i in 1:n^2
    # Plot a density heatmap on each subplot
    E_rounded = round(diag.values[i], sigdigits=4)
    plot!(p[i], localdensity(diag[i]), title="\$E_{$i} = $E_rounded\$", clims=clims, cbar=:none)
end

# The following 2 lines are kinda hacky; they draw one colorbar for all heatmaps
plot!(p[n^2+1], framestyle=:none)
scatter!([NaN], zcolor=[NaN], clims=clims, leg=:none, cbar=:right, background_subplot=:transparent, 
    framestyle=:none, inset=bbox(0.0, 0.05, 0.95, 0.9), subplot=n^2+2, c=:matter)
```

See more examples in the [Examples](@ref) section.