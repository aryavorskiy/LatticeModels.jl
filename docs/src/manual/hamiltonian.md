# Systems and Hamiltonians

This section describes how to define the Hamiltonian of a system. The system, in this context, is a lattice, a basis describing on-site degrees of freedom, and some additional parameters like temperature or chemical potential.

## Basic tight-binding Hamiltonian

First of all, we need to define the system. This is done by calling [`System`](@ref). For example, if we want to create non-interacting fermions with spin 1/2, zero chemical potential and temperature 1, we can do the following:

```@example 1
using LatticeModels, Plots
l = SquareLattice(4, 4);
sys = System(l, SpinBasis(1//2), μ = 0, T = 1)
H = tightbinding_hamiltonian(sys)
P = densitymatrix(H, info=false)
heatmap(localdensity(P))
```

In this case we have created a system with defined chemical potential and temperature. Note that you can use `mu` instead of `μ` if you prefer (or if you are using a non-UTF8 compatible editor). The [`tightbinding_hamiltonian`](@ref) function creates a tight-binding Hamiltonian for the given system, and the [`densitymatrix`](@ref) function calculates the density matrix.

The first two arguments of the `System` constructor are the lattice and the basis of the on-site degrees of freedom. 
The latter can be any `QuantumOptics.Basis` object - refer to the [QuantumOptics.jl documentation](https://qojulia.org/docs/) for more information. Omit the second argument if you don't have any on-site degrees of freedom.

Here are all parameters that can be passed to the `System` constructor:
- `mu` or `μ`: the chemical potential of the system.
- `N`: set this instead of `mu` to fix the number of particles - the chemical potential will be calculated automatically. Setting both `N` and `mu` will raise an error.
- `statistics`: the statistics of the particles, either `FermiDirac` or `BoseEinstein`. If not set, the statistics is `FermiDirac` if `mu` or `N` is set, and Gibbs otherwise (i.e. the system consist of one particle).
- `T`: the temperature of the system. The default is 0.

Note that you can pass the very same arguments directly to the `tightbinding_hamiltonian`. Also keyword arguments 
can be passed to the `densitymatrix` function - they will be used to evaluate the distribution function. These 
lines are equivalent to the previous example:

```julia
H = tightbinding_hamiltonian(l, SpinBasis(1//2))
P = densitymatrix(H, μ = 0, T = 1, info = false)
```

The `tightbinding_hamiltonian` function creates a basic tight-binding Hamiltonian defined by this formula:

```math
H = t_1 \sum_{\langle i, j \rangle} c_i^\dagger c_j + \text{h.c.} + 
    t_2 \sum_{\langle\langle i, j \rangle\rangle} c_i^\dagger c_j + \text{h.c.} + 
    t_3 \sum_{\langle\langle\langle i, j \rangle\rangle\rangle} c_i^\dagger c_j + \text{h.c.} +
    \ldots
```

where $c_i^\dagger$ and $c_i$ are the creation and annihilation operators on site $i$. The ``t_1``, ``t_2`` and ``t_3`` are responsible for the hopping between the nearest, next-nearest, etc. neighbors. Note that they can be factors or operators that act on the on-site degrees of freedom.

By default only the nearest-neighbor hopping is included, all other factors are zero. You can set them to any value you want using the `t1`, `t2`, `t3`, keyword arguments:

```@example 1
sz = sigmaz(SpinBasis(1//2))
H2 = tightbinding_hamiltonian(sys, t1 = sz, t2 = [1 0; 0 -1], t3 = -0.5)
P2 = densitymatrix(H2, info=false)
heatmap(localdensity(P2))
```

Here `sz` is the operator that projects the spin on the $z$ axis. Note that we can use a matrix instead, how we did with `t2`: in this case it will be interpreted as the hopping operator acting on the on-site degrees of freedom.

## The constructor function

## The operator builder

## Gauge fields