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

## Solvers

## TimeSequence