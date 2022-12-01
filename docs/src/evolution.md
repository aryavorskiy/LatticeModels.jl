We often  want to study the behavior of some quantum system in time-dependent conditions. We can use the unitary evolution operator to describe how the density matrix depends on time:

$$\mathcal{U}(t) = T\left\{ e^{\frac{1}{i\hbar} \int_{t_0}^t \hat{H}(\tau) d\tau} \right\},\hspace{0.5cm}
\mathcal{P}(t) = \mathcal{U}(t) \mathcal{P}_0 \mathcal{U}^\dagger (t)$$

## The evolution macro

This macro can be quite useful if your hamiltonian depends on time or if there are multiple hamiltonians in your experiment.
It avoids excessive computation in several cases automatically:
- If the hamiltonian does not change, the $\mathcal{U}(t, dt)$ operator will not be re-evaluated.
- If several wave functions or density matrices evolve using the same hamiltonian, neither the hamiltonian nor the evolution operator matrix will be re-evaluated.

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

The [`@evolution`](@ref) macro creates a scope where the hamiltonian and wave functions (or density matrices) are evaluated for the given time interval. 
It takes two arguments: a braced list with evolution specifiers and a for-loop that iterates over the time interval:

```@example env
l = SquareLattice(10, 10)
τ = 30
a = Animation()
B = 0.01
P0 = filled_projector(spectrum(Chern(l, 0)))

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

## Compatibility

This macro is designed to be compatible with arbitrary array types. This means that hamiltonian expressions or initial states needn't to be `LatticeArray`s, but they can be of any array type instead, as long as it supports linear algebra operations such as matrix multiplication.

Here are the functions that must be defined for the desired array type:

- Equality operator: `==`
- Basic arithmetic functions: `+`, `-`, `*`, `adjoint`
- The identity matrix: `one(A)`
- The matrix exponent `exp(A)`
    - If it is not possible to implement this function you can set the `k` keyword argument (see below) to calculate the matrix exponent as a partial sum of Taylor series.

!!! warning
    If you use `LatticeArray`s, you still have to make sure these functions and operators are defined for the underlying array type.

## Performance

The evolution macro avoids calculations where possible to improve performance. It is important to know how it achieves this result:

- If the hamiltonian matrix has not changed between two iterations and the time step remained approximately the same, $\mathcal{U}(t, dt)$ will not be re-evaluated, because the matrix exponent is the most time-consuming operation compared to others like matrix multiplication or addition.
  - The relative tolerance used to check if the time step has changed can be set via `rtol` keyword.
- If the hamiltonian expression does not explicitly depend on the loop variable (`t` in the example), it will be considered constant and evaluated only once at the beginning. Otherwise it will be evaluated on every iteration in the loop scope.
- If several states evolve according to the same hamiltonian, both the hamiltonian and the $\mathcal{U}(t, dt)$ evolution operator will be evaluated only once per iteration.

To improve performance with small time intervals you can calculate the matrix exponent as a Taylor polynomial. Its order can be set via `k` keyword.

By default the macro shows a progress informer that shows the task progress, the estimated time remaining and the fraction of time that
was spent to perform evolution. To disable it add `show_progress=false` to the keyword arguments.

Keyword assignments should be placed before the rules list:

```julia
@evolution k=2 rtol=1e-6 show_progress=false {...} for t in ...
```