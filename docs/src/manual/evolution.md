# Evolution

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

We got a pertty plot of a Gaussian wavepacket moving to the right. Let us discuss what happened here in-depth.

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

Note that we used a different notation to access the state and the time. You can use any of these two, whichever you find more convenient in your situation:
- Fields: `psi = moment.state; Ht = moment.H; t = moment.t`
- Unpacking: `psi, Ht, t = moment`

## Time-dependent Hamiltonian

If the Hamiltonian is time-dependent, you can pass a function that returns the Hamiltonian at the given time. This function should accept a single argument — the time. An example of a time-dependent Hamiltonian is the one with a magnetic field that changes in time from the [Examples section](@ref "examples"):

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

    Note that you can also use [`QuantumOptics.TimeDependentOperator`](https://docs.qojulia.org/timeevolution/timedependent-problems/#Time-dependent-operators) interface — this might be more convenient, as it eliminates extra allocations, and exponent caching for them is disabled by default.

## Solvers

## TimeSequence