## Usage

An [`AbstractCurrents`](@ref) object is a lazy object that can calculate any current-like value between any pair of sites. 
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
    field := FluxField(1, (5.5, 5.5))
end
P = filled_projector(spectrum(H), -0.5)

curr = DensityCurrents(H, P)                    # Create Currents object
heatmap(site_density(P))
plot!(curr, arrows_scale=25, arrows_rtol=0.1, color=:blue)    # Quiver-plot the currents
```

What happened here? The formula for the density current from site $i$ to site $j$ is $J_{ij} = \text{tr}(-i \hat{h}_{ij} \hat{c}^\dagger_i \hat{c}_j \hat{\rho} + h. c.) = 2 \text{Im tr}(\hat{h}_{ij} \hat{c}^\dagger_i \hat{c}_j \hat{\rho})$. 
The `curr` object contains this formula inside and when this object is passed to the `plot` function, all needed currents are evaluated.

A current between each pair of sites is shown as an arrow directed from one to another with its length proportional to the strength of the current. The `arrows_scale` keyword argument scales all arrows by given factor, while `arrows_rtol` hides all arrows that are shorter than some fraction of the distance between the sites.

You can display currents between certain sites by using a boolean-typed `LatticeValue` and bracket-notation:

```@example env
x, y = coord_values(l)
sub_curr = curr[x .< y]
plot!(sub_curr, arrows_scale=25, arrows_rtol=0.1, color=:green)
```

Here all currents **between** sites in the upper-left coordinate triangle are marked green.

## Interface

It is quite likely that you might want to define your own type of currents. All you need to do is inherit the `AbstractCurrents` type and define two functions:

- `lattice(::MyCurrents)` must return the lattice which the currents are defined on
- `Base.getindex(::MyCurrents, i::Int, j::Int)` must return the current between sites with indices `i` and `j`. Note that the function must be skew-symmetric, e. g. `curr[i, j] == -curr[j, i]`.

## Materialized currents

An `AbstractCurrents` is a lazy object. This allows to avoid excessive computation of site-to-site currents, but the computations that are needed will be repeated every time when we access that object; also abstract currents cannot be normally stored into a `TimeSequence` (more precisely, you won't be able to differentiate or integrate such currents over time). That's where the `MaterializedCurrents` come in, having all their values stored explicitly in an array.

To convert any type of currents to `MaterializedCurrents`, simply use the [`materialize`](@ref) function. You can avoid evaluating some currents (for example, if you know beforehand that they must be zero) by passing a lambda as a first argument (or with `do`-syntax): it must take the `Lattice` and two `LatticeSite`s and return whether the current between these sites must be evaluated.

You can find it similar to the selector function we used back in [Hopping operators](@ref Hopping-operators), which indeed is.
You may find the following selector functions useful:

- Passing a `PairSet` produced by the [`bonds`](@ref) function will keep only the currents between adjacent sites.
- [`pairs_by_distance`](@ref) will allow you to select pairs of sites depending on the distance between them. 

## Mapping currents

In some cases we want to find out how currents depend on some lattice properties: for example, the distance between sites.
In such case, the [`map_currents`](@ref) function can be quite helpful.

Let's find the mean and the standard deviation for currents between sites given the distance between them:

```@example env
using LinearAlgebra, Statistics

dist, adcurr = map_currents(
    curr, 
    reduce_fn=(x -> [mean(abs.(x)), std(abs.(x))]),
    sort=true
) do l, site1, site2
    norm(site1.coords - site2.coords)
end

acurr, dcurr = eachcol(adcurr)
scatter(dist, acurr, err=dcurr, xlims=(0, 14))
```

What happened here? The `map_currents` function found the distance and the current between each pair of sites. Then for each distance between sites it found the mean and standard deviation for the absolute value of the currents in such pairs, and stored it column-wise in a matrix automatically.
In the next line we extracted the mean and standard deviation into separate lists, and plotted the obtained data.

From this picture we can see that there are no density currents between non-adjacent sites, as one must have expected.