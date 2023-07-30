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
- Smart unitary evolution macro reducing excessive computations where possible.
- Supports visualization via [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- Compatible with [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Usage example

This simple code plots density heatmaps for lowest energy states of a square tight-binding lattice.

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
    plot!(p[i], site_density(diag[i]), title="\$E_{$i} = $E_rounded\$", clims=clims, cbar=:none)
end

# The following 2 lines are kinda hacky; they draw one colorbar for all heatmaps
plot!(p[n^2+1], framestyle=:none)
scatter!([NaN], zcolor=[NaN], clims=clims, leg=:none, cbar=:right, background_subplot=:transparent, 
    framestyle=:none, inset=bbox(0.0, 0.05, 0.95, 0.9), subplot=n^2+2)
```

See more examples in the [Tutorial](@ref) section.