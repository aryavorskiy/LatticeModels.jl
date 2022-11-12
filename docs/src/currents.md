## Usage

A [`AbstractCurrents`](@ref) object is a lazy object that can calculate any current-like value between any pair of sites. 
Usage is pretty simple:

```@setup env
using LatticeModels, Plots
```

```@example env
l = SquareLattice(10, 10)
H = @hamiltonian begin   
    lattice := l
    @hop axis = 1
    @hop axis = 2
    field := FluxField(1, (0, 0))
end
P = filled_projector(spectrum(H), -0.5)

curr = DensityCurrents(H, P)                    # Create Currents object
heatmap(ptrace(P) .|> real)
plot!(curr, arrows_scale=25, arrows_rtol=0.1, color=:blue)    # Quiver-plot the currents
```

What happened here? The formula for the density current from site $i$ to site $j$ is $J_{ij} = \text{tr}(-i \hat{h}_{ij} \hat{c}^\dagger_i \hat{c}_j \hat{\rho} + h. c.) = 2 \text{Im tr}(\hat{h}_{ij} \hat{c}^\dagger_i \hat{c}_j \hat{\rho})$. 
The `curr` object contains this formula inside - and when this object is passed to the `plot` function, all needed currents are evaluated.

A current between each pair of sites is shown as an arrow directed from one to another with its length proportional to the strength of the current. The `arrows_scale` keyword argument scales all arrows by given factor, while `arrows_rtol` hides all arrows that are shorter than some fraction of the distance between the sites.

You can display currents between certain sites by using a boolean-typed `LatticeValue` and bracket-notation:

```@example env
x, y = coord_values(l)
sub_curr = curr[x .< y]
plot!(sub_curr, arrows_scale=25, arrows_rtol=0.1, color=:green)
```

## Interface

It is quite likely that you might want to define your own type of currents. All you need to do is inherit the `AbstractCurrents` type and define two functions:

- `lattice(::MyCurrents)` must return the lattice which the currents are defined on
- `currents_lambda(::MyCurrents)` must return a lambda which takes two integer site indices and returns the current between these two sites. Note that the function must be skew-symmetric, e. g. `curr_lambda(i, j) == -curr_lambda(j, i)`

## Materialized currents

An `AbstractCurrents` is a lazy object - this allows to avoid excessive computation, but the computations that are needed will be repeated every time when we use the object. That's where the `MaterializedCurrents` come in, having all their values stored explicitly in an array.

To convert any type of currents to `MaterializedCurrents`, simply use the [`materialize`](@ref) function. You can avoid evaluating some currents (for example, if you know beforehand that they must be zero) by passing a lambda as a first argument (or with `do`-syntax): it must take the `Lattice` and two integer indices and return whether the current between these sites must be evaluated.

You can find it similar to the selector function we used back in [Hopping operators](@ref Hopping-operators), which indeed is.
You may find the following selector functions useful:

- [`pairs_by_adjacent`](@ref) will keep only the currents between adjacent sites.
- [`pairs_by_distance`](@ref) will allow you to select pairs of sites depending on the distance between them. 