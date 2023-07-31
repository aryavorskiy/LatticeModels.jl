We often  want to study the behavior of some quantum system in time-dependent conditions. We can use the unitary evolution operator to describe how the density matrix depends on time:

$$\mathcal{U}(t) = T\left\{ e^{\frac{1}{i\hbar} \int_{t_0}^t \hat{H}(\tau) d\tau} \right\},\hspace{0.5cm}
\mathcal{P}(t) = \mathcal{U}(t) \mathcal{P}_0 \mathcal{U}^\dagger (t)$$

## The evolution macro

This macro can be quite useful if your hamiltonian depends on time or if there are multiple hamiltonians in your experiment.
It avoids excessive computation in several cases automatically:
- If the hamiltonian does not change, the $\mathcal{U}(t, dt)$ operator will not be re-evaluated.
- If several wavefunctions or density matrices evolve using the same hamiltonian, neither the hamiltonian nor the evolution operator matrix will be re-evaluated.

Let us define a function that generates a Chern insulator hamiltonian:

```@example env
using LatticeModels, Plots
σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
Chern(l, B) = @hamiltonian begin
    lattice := l
    field := LandauField(B)
    dims_internal := 2
    @diag σ[3]
    @hop axis=1 (σ[3] - im * σ[1]) / 2
    @hop axis=2 (σ[3] - im * σ[2]) / 2
end
nothing # hide
```

The [`@evolution`](@ref) macro creates a scope where the hamiltonian and wavefunctions (or density matrices) are evaluated for the given time interval. 
It takes two arguments: a braced list with evolution specifiers and a for-loop that iterates over the time interval:

```@example env
l = SquareLattice(10, 10)
τ = 30
a = Animation()
B = 0.01
P0 = densitymatrix(diagonalize(qwz(l, 0)))

@evolution {
    H := Chern(l, B * min(t, τ) / τ),
    P0 --> H --> P
} for t in 0:0.1:2τ
    cur = DensityCurrents(H, P)
    heatmap(site_density(P), title="Local density, t = $t", clims=(0.9, 1.1))
    plot!(cur, arrows_scale=20, color=:blue)
    frame(a)
end
gif(a, "animation.gif")
```

Let us make it clear what an evolution specifier is. In fact, there are two possible variants:

- It is an *evolution rule* describing the initial state, the hamiltonian it will evolve accroding to, and the name of the variable to write the result. These three arguments must be written in a chain delimited by `-->` like this: 
  
  `P0 --> Chern(l, B * min(t, τ) / τ) --> P`

  Note that you can use a hamiltonian alias instead of an expression if it was previously defined.
- It is an *alias declaration*, which means that a certain expression will be evaluated on every iteration and assigned to a variable with a given name. For example, `H := Chern(l, B * min(t, τ) / τ)` evaluates the Chern insulator hamiltonian and writes the result to the `H` variable.

### Compatibility

This macro is designed to be compatible with arbitrary array types. This means that hamiltonian expressions or initial states needn't to be `LatticeArray`s, but they can be of any array type instead, as long as it supports linear algebra operations such as matrix multiplication.

Here are the functions that must be defined for the desired array type:

- Equality operator: `==`
- Basic arithmetic functions: `+`, `-`, `*`, `adjoint`
- The identity matrix `one(A)`
- The matrix exponent `exp(A)`
    - If it is not possible to implement this function you can set the `k` keyword argument (see below) to calculate the matrix exponent as a partial sum of Taylor series.
- If the inverse matrix `inv(A)` function is defined, you can use the Padé approximant for the matrix exponent, which is generally more accurate.

!!! warning
    If you use `LatticeArray`s, you still have to make sure these functions and operators are defined for the underlying array type.

### Performance

The evolution macro avoids calculations where possible to improve performance. It is important to know how it achieves this result:

- If the hamiltonian matrix has not changed between two iterations and the time step remained approximately the same, $\mathcal{U}(t, dt)$ will not be re-evaluated, because the matrix exponent is the most time-consuming operation compared to others like matrix multiplication or addition.
  - The relative tolerance used to check if the time step has changed can be set via `rtol` keyword.
- If the hamiltonian expression does not explicitly depend on the loop variable (`t` in the example), it will be considered constant and evaluated only once at the beginning. Otherwise it will be evaluated on every iteration in the loop scope.
- If several states evolve according to the same hamiltonian, both the hamiltonian and the $\mathcal{U}(t, dt)$ evolution operator will be evaluated only once per iteration.

To improve performance with small time intervals you can calculate the matrix exponent as a Taylor polynomial. Its order can be set via `k` keyword.

For more precise calculation add `pade=true` - this will enable matrix exponent calculation via Padé approximant. Note that this formula requires finding an inverse matrix, so this option is not compatible with sparse matrices.

By default the macro shows a progress informer that shows the task progress, the estimated time remaining and the fraction of time that
was spent to perform evolution. To disable it add `show_progress=false` to the keyword arguments.

Keyword assignments should be placed before the rules list:

```julia
@evolution k=2 rtol=1e-6 show_progress=false {...} for t in ...
```

## Time sequences

A `TimeSequence` is a struct that stores information about how some value changed during time. It simplifies working with time-dependent values by allowing you to run the computation pass only once and re-evaluate all processing and visualization code as much as you want.

Here is an example:

```@example env
P0 = densitymatrix(diagonalize(qwz(l)))
density_rec = TimeSequence()
deriv_rec = TimeSequence()

@evolution {
    H := Chern(l, B * min(t, τ) / τ),
    P0 --> H --> P
} for t in 0:0.1:2τ
    insert!(density_rec, t, site_density(P))
    insert!(deriv_rec, t, site_density(-im * (H * P - P * H)))
end

site = l[50]
p = plot(layout=(2,1))
plot!(p[1], density_rec[site], lab="p")

# Compare computed time derivative with Heisenberg equation
plot!(p[2], differentiate(density_rec[site]), lw=5, alpha=0.3, lab="dp/dt")
plot!(p[2], deriv_rec[site], lab="Heisenberg")
```

### Indexing and iteration

`TimeSequence`s support two kinds of indexing:

- Calling the record with a numeric value selects the time and returns the object in nearest snapshot. 
  Calling it with two numeric values yields a new `TimeSequence` with all timestamps between given values.
- Indexing it with square brackets will apply this index to all snapshots. The output will also be a `TimeSequence`.

Let's see how it works:

```@example env
p = plot(layout=(1,3), size=(800, 250))
plot!(p[1], density_rec(10))                        # The snapshot with time nearest to 10
typeof(density_rec(10))                             # LatticeValue{Float64}
plot!(p[2], density_rec[l[25]])                     # The value on the 25-th site depending on time 
typeof(density_rec[l[25]])                          # Dict{Float64, Float64}
plot!(p[3], (density_rec(15, 25) |> diff)[l[25]])   # Derivative of the value by time where 15 ≤ t ≤ 25
```

Note that `TimeSequence` impements the `AbstractDict` and therefore is iterable. The iterator yields `timestamp => value` pairs:

```julia
for (t, P) in density_rec
    # Do what you need with t, P
end
```

### Differentiation and integration

You can also automatically differentiate and integrate the `TimeSequence` with [`differentiate`](@ref) and [`integrate`](@ref) functions (or their in-place alternatives [`differentiate!`](@ref) and [`integrate!`](@ref)). However, the type of the stored values must implement the following methods:
- `Base.copy`: the non in-place methods simply call their in-place versions on a copy of the original TimeSequence.
- `Base.zero`: this function will be called before integrating.
- `LinearAlgebra.axpby!` to implement basic in-place math.

Still, types in this library (`LatticeValue`, `LatticeArray`, `MaterializedCurrents`), numbers and arrays are compatible with differentiation and integration out of the box.

!!! warning
    Integrating or differentiating may produce unexpected results or crash if applied to a `TimeSequence` with elements of different type or defined on a different [`Lattice`](@ref)/[`Basis`](@ref). **Lattice or basis match is not checked when constructing a `TimeSequence`!**