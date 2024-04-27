# [Evolution](@id evolution_chapter)

This chapter describes several ways to work with time evolution of quantum systems on lattices. It features [`Evolution`](@ref) — a powerful struct that represents the time evolution of a quantum system according to the Schrödinger equation, and `TimeSequence`, which is used to store and process time-dependent data.

## Basics

The `Evolution` struct is used to represent the time evolution of a quantum system. It contains the Hamiltonian and the initial state (or states) of the system. To evaluate the time evolution, you can write this:

```@example 1
using LatticeModels, Plots
lat = Chain(200)

# Define the initial state: a Gaussian wavepacket
x = coordvalue(lat, :x)
psi0 = ketstate(@. exp(-0.02 * (x - 50) ^ 2 - im * 0.3x))

# Define the Hamiltonian
H = tightbinding_hamiltonian(lat)

# Create the plot
plot(title="Gaussian wavepacket evolution", xlabel="x", yticks=:none)
# Empty plots for the legend
plot!([NaN], c = :blue, label="Re") 
plot!([NaN], c = :red, label="Im")  

ev = Evolution(H, psi0)
for moment in ev(0:5:50)
    psi, Ht, t = moment # Unpack to get the state, Hamiltonian and time
    plot!(real.(psi.data) .+ t / 3, c=:blue, lab="")
    plot!(imag.(psi.data) .+ t / 3, c=:red, lab="")
end
plot!()
```

We got a pretty plot of a Gaussian wavepacket moving to the right. Let us discuss what happened here in-depth.

We set the Hamiltonian to be a constant operator. Since we perform the evolution in even time steps, and the Hamiltonian is not time-dependent, we can precompute the exponential of the Hamiltonian ``\hat{U} = e^{-i H t}`` and use it in the evolution. This is done automatically here.

The `Evolution` struct is not an iterator. Instead, it is a callable object that returns an iterator. This is done to allow for more flexibility, as you can see below.

The `Evolution` struct always yields stateful iterators. This means that the `psi` in the loop is actually the same object, and it is updated on each iteration. This is done to reduce memory allocations and improve performance. Do not edit the states in the loop, as it will affect the results.

This allows continuing the evolution from the last state. Let's stick to the previous example and add some more iterations:

```@example 1
for moment in ev(50:5:100)
    t = moment.t        # Get the time
    psi = moment.state  # Get the state
    # moment.H is the Hamiltonian, but it is not used here
    plot!(real.(psi.data) .+ t / 3, c=:blue, lab="")
    plot!(imag.(psi.data) .+ t / 3, c=:red, lab="")
end
plot!()
```

!!! tip
    If you do not need to iterate over a single `Evolution` object multiple times, you can use a shorthand. This notation:

    ```julia
    for moment in Evolution(H, psi0, timedomain=0:5:50)
        # ...
    end
    ```
    
    is equivalent to this:

    ```julia
    ev = Evolution(H, psi0)
    for moment in ev(0:5:50)
        # ...
    end
    ```

Note that we used a different notation to access the state and the time. You can use any of these two, whichever you find more convenient in your situation:
- Fields: `psi = moment.state; Ht = moment.H; t = moment.t`
- Unpacking: `psi, Ht, t = moment`

## Multiple initial states

The `Evolution` struct can also accept multiple initial states. Let's create two Gaussian wavepackets and evolve them:

```@example 1
# Define the initial states: two Gaussian wavepackets
psi0_1 = ketstate(@. exp(-0.02 * (x - 50) ^ 2 - im * 0.3x))
psi0_2 = ketstate(@. exp(-0.02 * (x - 150) ^ 2 + im * 0.3x))

p = plot(title="Two Gaussian wavepackets evolution", xlabel="x", yticks=:none)

ev = Evolution(H, psi0_1, psi0_2)     # Pass the initial states as keyword arguments
for moment in ev(0:5:50)
    psi1 = moment[1]
    psi2 = moment[2]
    t = moment.t
    plot!(abs2.(psi1.data) .+ t / 3, c=:black, lab="")
    plot!(abs2.(psi2.data) .+ t / 3, c=:green, lab="")
end
plot!()
```

Note that accessing `moment.state` will throw an error here, because there is more than one state. Instead, you can access all them as a tuple with `moment.states` or each one separately with `moment[1]`, `moment[2]`, etc. Also you can still unpack the `moment` as `psi1, psi2, Ht, t = moment`, which will do the very same thing.

The `Evolution` struct also allows assigning aliases to the states. You can do this by passing the initial states as keyword arguments. This code will produce the same result as in the previous example:

```julia
p = plot(title="Two Gaussian wavepackets evolution", xlabel="x", yticks=:none)

ev = Evolution(H, psi1=psi0_1, psi2=psi0_2)
for moment in ev(0:5:50)
    t = moment.t
    plot!(abs2.(moment.psi1.data) .+ t / 3, c=:black, lab="")
    plot!(abs2.(moment.psi2.data) .+ t / 3, c=:green, lab="")
end
plot!()
```

`moment.states` will be a named tuple with the states in this case. This way you can access the states as `moment.psi1`, `moment.psi2`, etc. You can still use indices `moment[1]`, `moment[2]` or unpack the `moment` as `psi1, psi2, Ht, t = moment`. The choice is yours.

## Time-dependent Hamiltonian

If the Hamiltonian is time-dependent, you can pass a function that returns the Hamiltonian at the given time. This function should accept a single argument — the time. An example of a time-dependent Hamiltonian is the one with a magnetic field that changes in time from the [Examples section](@ref "Examples"):

```@example 2
using LatticeModels, Plots
l = TriangularLattice(Circle(10), !Circle(5))
removedangling!(l)
h(B) = tightbinding_hamiltonian(l, field=PointFlux(B))
diag = diagonalize(h(0))

P_0 = densitymatrix(diag, mu = 0)
τ = 10
ht(t) = h(0.1 * min(t, τ) / τ)
ev = Evolution(ht, P_0)

anim = @animate for moment in ev(0:0.1:2τ)
    ρ, Ht, t = moment
    curr = DensityCurrents(Ht, ρ)
    plot(curr, title="t = $t", clims=(0, 0.005), size=(1000, 700))
end
gif(anim, "currents.gif")
```

Here we applied magnetic field that changes in time. The flux is increased linearly from 0 to 0.1 in 10 time units, and then stays constant. The animation shows the density currents in the lattice at each time step.

To make the Hamiltonian time-dependent in the previous example, we defined a function `ht(t) = h(0.1 * min(t, τ) / τ)`. This function returns the Hamiltonian with the flux `0.1 * min(t, τ) / τ` at the given time `t`. The `Evolution` struct then uses this function to calculate the Hamiltonian at each time step.

The exponent for the time-dependent Hamiltonian is cached for each time step. This means that on the second half of the evolution interval, the exponent is calculated only once, and then reused for each time step. This improves performance. To improve it even more, you can pre-calculate the static Hamiltonian:

```julia
const H1 = h(0.1)   # Make it constant to enforce type stability
ht(t) = t > τ ? H1 : h(0.1 * t / τ)
```

!!! note
    The `Evolution` checks the Hamiltonian equality by storing a reference to the last Hamiltonian. Thes means that if you return the same object each time in the Hamiltonian function, the exponent will not be recalculated, even if you changed the Hamiltonian. To force the recalculation, wrap the return value into a `Ref`.

    Note that you can also use [`QuantumOptics.TimeDependentOperator`](https://docs.qojulia.org/timeevolution/timedependent-problems/#Time-dependent-operators) interface — this might be more convenient, as it eliminates extra allocations more efficiently, and exponent caching for them is disabled by default.

## Solvers

The `Evolution` struct can use several solvers to calculate the time evolution. You can set the solver to use explicitly this way:

```julia
ev = Evolution(CachedExp(threshold=1e-12, nztol=1e-15), H, psi0)
```

The [`CachedExp`](@ref) solver is the default option, and it works in a very simple way: at each time step it calculates the exponent of the Hamiltonian ``\hat{U} = e^{-i H \Delta t}`` and applies it to the state. It also caches the exponent to avoid recalculating it if the Hamiltonian and the timestep are the same.

The exponent is calculated using a custom-implemented scaling and squaring algorithm. Contrary to the `exp` function from the `LinearAlgebra` package, it works with sparse matrices and GPU arrays also, while performance on dense matrices is comparable. However, some precision loss is possible, especially for large time steps, so you can adjust the tolerance for your case with these parameters:

- `threshold` — the threshold for the error in the exponent calculation. Limits the number of iterations in the Taylor series expansion. The default value is `1e-10`.
- `nztol` — the tolerance for the zero elements in the Hamiltonian. Works only for sparse matrices — the elements with the absolute value less than `nztol` are considered zero. The default value is `1e-14`.

The `CachedExp` solver is universal and works well in most cases. However, it is not the most efficient one. In some specific cases you might want to use a more specialized solver. For example, the [`KrylovKitExp`](@ref) solver can be faster if the Hamiltonian is large, sparse and dependent on time. It uses the `exponentiate` function from the `KrylovKit` package to calculate the exponent.

Let's run some performance tests:

```@example 3
using LatticeModels
l = SquareLattice(40, 40)
h(B) = tightbinding_hamiltonian(l, field=LandauGauge(B))

# The initial state is the ground state of the Hamiltonian with B = 0
psi = groundstate(h(0))
τ = 10
ht(t) = h(0.1 * min(t, τ) / τ)
ts = 0:0.1:2τ
site = l[!, x = 20, y = 20]

println("CachedExp:")
# We can evaluate the density at x = 20, y = 20 in a one-liner
@time dens1 = [localdensity(moment.state)[site] for moment in 
    Evolution(CachedExp(), ht, psi, timedomain=ts, showprogress=false)]
@time dens1 = [localdensity(moment.state)[site] for moment in 
    Evolution(CachedExp(), ht, psi, timedomain=ts, showprogress=false)]
@time dens1 = [localdensity(moment.state)[site] for moment in 
    Evolution(CachedExp(), ht, psi, timedomain=ts, showprogress=false)]

println("\nKrylovKitExp:")
@time dens2 = [localdensity(moment.state)[site] for moment in 
    Evolution(KrylovKitExp(), ht, psi, timedomain=ts, showprogress=false)]
@time dens2 = [localdensity(moment.state)[site] for moment in 
    Evolution(KrylovKitExp(), ht, psi, timedomain=ts, showprogress=false)]
@time dens2 = [localdensity(moment.state)[site] for moment in 
    Evolution(KrylovKitExp(), ht, psi, timedomain=ts, showprogress=false)]
nothing # hide
```

The first attempt took some time because of the precompilation required, but overall the `KrylovKitExp` solver turned out to be faster in this case. Let's plot the results to make sure they are the same:

```@example 3
using Plots
plot(ts, dens1, label="CachedExp", xlabel="t", ylabel="ρ(x=20, y=20)", 
    ylims=(0, NaN), size=(800, 400))
plot!(ts, dens2, label="KrylovKitExp")
```

A significant drawback of `KrylovKitExp` solver is that it works only with wavefunctions — it does not support density matrices.

!!! note
    Generally, the `Evolution` struct supports only Schrödinger/von Neumann equations. If you need to solve the master equation, you can use [the functions from the `QuantumOptics` package](https://docs.qojulia.org/timeevolution/master/).

    However, a solver for the master equation can be implemented on top of the `Evolution` struct using this `EvolutionSolver` interface. Please file an issue if you are interested in this feature.

## TimeSequence

The [`TimeSequence`](@ref) struct is used to store and manipulate time-dependent data. It is a dictionary-like object that maps time points to values. The values can be any type, including arrays, matrices, and even lattice-specific types like `LatticeValue`.

Let's calculate the evolution of a ground state of a tight-binding model after a magnetic field is adiaiabatically turned on. We will store the local density at each time step and use it to plot the local density depending on time, as well as its time derivative and integral over time:

```@example 4
using LatticeModels, Plots
l = SquareLattice(20, 20)
h(B) = tightbinding_hamiltonian(l, field=LandauGauge(B))

# The initial state is the ground state of the Hamiltonian with B = 0
psi = groundstate(h(0))
τ = 10
ht(t) = h(0.1 * min(t, τ) / τ)
densities = TimeSequence{LatticeValue}()

for moment in Evolution(ht, psi, timedomain=0:0.1:2τ)
    densities[moment.t] = localdensity(moment.state)
end

plot(densities[!, x = 10, y = 10], label="ρ(t) (bulk)")
plot!(densities[!, x = 1, y = 10], label="ρ(t) (edge)")
```

!!! tip
    The very same result could be achieved without writing the loop, using the `TimeSequence` constructor with a function. This way it will automatically iterate over the `Evolution` object and store the results:

    ```julia
    densities = TimeSequence(Evolution(ht, psi), 0:0.1:2τ) do moment
        localdensity(moment.state)
    end
    ```

Note that we can index a `TimeSequence` of `LatticeValue`s the same way as a `LatticeValue` itself, which will produce another `TimeSequence` with the same time domain. This is because `TimeSequence` forwards almost all the indexing operations to the values it stores.

This type also supports time indexing and slicing. You can combine these types of indexing by using the `t` keyword for the time axis — all other arguments will be treated as spatial coordinates. Let's plot the local density at the center of the lattice and the local density at the edge of the lattice at the same time point:

```@example 4
p = plot(layout=2, size=(1000, 500))
plot!(p[1], densities[0.5])
plot!(p[2], densities[!, x=10, y=10, t=0.5 .. 3.5], label="ρ(t) (bulk)")
plot!(p[2], densities[!, x=1, y=10, t=5 .. 10], label="ρ(t) (edge)")
```

Another very useful feature of the `TimeSequence` is that it can be integrated and differentiated over time using the [`integrate`](@ref) and [`differentiate`](@ref) functions. Let's plot the time derivative and integral of the local density in the bulk of the lattice:

```@example 4
densities_bulk = densities[!, x = 10, y = 10]
plot(densities_bulk, label="ρ(t) (bulk)")
plot!(differentiate(densities_bulk), label="dρ(t)/dt (bulk)")
plot!(integrate(densities_bulk), label="∫ρ(t)dt (bulk)")
```