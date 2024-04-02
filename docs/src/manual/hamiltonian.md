# Systems and Hamiltonians

This section describes how to define the Hamiltonian of a system. The system, in this context, is a lattice, a basis describing on-site degrees of freedom, and some additional parameters like temperature or chemical potential.

## Defining the system

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
The latter can be any `QuantumOptics.Basis` object — refer to the [QuantumOptics.jl documentation](https://docs.qojulia.org/quantumobjects/bases/) for more information. Omit the second argument if you don't have any on-site degrees of freedom.

Here are all parameters that can be passed to the `System` constructor:
- `mu` or `μ`: the chemical potential of the system.
- `N`: set this instead of `mu` to fix the number of particles — the chemical potential will be calculated automatically. Setting both `N` and `mu` will raise an error.
- `statistics`: the statistics of the particles, either `FermiDirac` or `BoseEinstein`. If not set, the statistics is `FermiDirac` if `mu` or `N` is set, and Gibbs otherwise (i.e. the system consist of one particle).
- `T`: the temperature of the system. The default is 0.

Note that you can pass the very same arguments directly to the `tightbinding_hamiltonian`. Also keyword arguments 
can be passed to the `densitymatrix` function — they will be used to evaluate the distribution function. These 
lines are equivalent to the previous example:

```julia
H = tightbinding_hamiltonian(l, SpinBasis(1//2))
P = densitymatrix(H, μ = 0, T = 1, info = false)
```

## Basic tight-binding Hamiltonian

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

There are many cases where more complex Hamiltonians are needed. In this case, you can use the [`construct_hamiltonian`](@ref) function. This function takes a system and several arguments, each describing a term in the Hamiltonian. As an example, let's consider the same Hamiltonian as in the previous example:

```@example 2
using LatticeModels, Plots
l = SquareLattice(4, 4)
spin = SpinBasis(1//2)
sys = System(l, spin, μ = 0, T = 1)

sz = sigmaz(spin)
H = construct_hamiltonian(sys, 
    sz => NearestNeighbor(1),
    [1 0; 0 -1] => NearestNeighbor(2),
    -0.5 => NearestNeighbor(3))
P = densitymatrix(H, info=false)
heatmap(localdensity(P))
```

Basically, every term is a pair of an operator acting on the on-site degrees of freedom and an object describing the lattice part. The lattice part can be one of the following:

- **On-site term**: ``\sum_i t_i c_i^\dagger c_i``. In this case, the operator is applied to every site.
    - A single `site` describes a single ``t_i c_i^\dagger c_i`` term. A convenient way to create an on-site potential.
    - A `LatticeValue` object describes the `t_i` coefficients for every site.
- **Hopping term**: ``\sum_{\langle i, j \rangle} t_{ij} c_i^\dagger c_j + h.c.``. Terms like this are used to describe the hopping between different sites.
    - A single `site1 => site2` pair describes a single hopping term between `site1` and `site2`.
    - An `AbstractBonds` object describes the hopping topology. See the [Adjacency and boundary conditions](@ref) section for more information on them.
- **Arbitrary operator**: any `QuantumOptics.Operator` object can be passed as a term, in which case it will be automatically embedded into the Hamiltonian.
    - Note that a matrix can be passed as an operator, in which case it will be interpreted as one acting on the on-site degrees of freedom. To act on the lattice, use `Operator(l, matrix)` — this will create a lattice operator that acts on the lattice `l` with the given matrix.

To showcase the flexibility of the `construct_hamiltonian` function, let's consider a more complex example. We will create a Hamiltonian of the QWZ model of a topological insulator (on a square lattice).

The QWZ model Hamiltonian is defined by the following formula:

```math
\hat{H} =
\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i +
\sum_i^\text{sites} \left(
c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i +
c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i +
h. c. \right)
```

We will take $m_i = 1$ and add some disorder to the on-site term. Also we will add electric field in the $x$ direction, local on-site potential on site at $(2, 2)$ and a new hopping term between the next-nearest neighbors.

To do this, we need to create the Hamiltonian like this:

```@example 3
using LatticeModels, Plots
l = SquareLattice(4, 4)
spin = SpinBasis(1//2)
sys = System(l, spin, μ = 0, T = 1)

sx, sy, sz = sigmax(spin), sigmay(spin), sigmaz(spin)
site = l[!, x = 2, y = 2]
x = LatticeValue(l, :x)
E = 0.1                                     # Electric field strength

H = construct_hamiltonian(sys,
    sz => 1 .+ 0.1 .* randn(l),             # Random on-site potential
    (sz - im * sx) / 2 => Bravais[1, 0],    # x-direction hoppings
    (sz - im * sy) / 2 => Bravais[0, 1],    # y-direction hoppings
    0.5 => site,                            # Local potential on site (2, 2)
    E * x,                                  # Electric field in the x direction
    0.1 => NearestNeighbor(2),              # Next-nearest neighbor hopping
    )

P = densitymatrix(H, info=false)
heatmap(localdensity(P))
```

Note that the [`qwz`](@ref) function can be used to create the QWZ Hamiltonian in a more convenient way. The `construct_hamiltonian` function is more flexible and can be used to create any Hamiltonian you want.

!!! note
    The `construct_hamiltonian` function creates a `Hamiltonian` object, which is a `QuantumOptics.Operator` object,
    but with some additional info like particle statistics or chemical potential. You can use it as a regular 
    operator, but if you encounter any problems, you can always convert it to a regular operator using the `Operator`
    function or by calling `construct_operator` instead of `construct_hamiltonian`.

## The operator builder

The `construct_hamiltonian` function is a convenient way to create Hamiltonians, but it is sometimes more convenient to set all the individual ``c_i^\dagger c_j`` terms manually. This can be done using the `OperatorBuilder` object.

Here is an example of how to create the same Hamiltonian as in the previous example using the `OperatorBuilder`:

```@example 3
using LatticeModels, Plots
l = SquareLattice(4, 4)
spin = SpinBasis(1//2)
sys = System(l, spin, μ = 0, T = 1)

sx, sy, sz = sigmax(spin), sigmay(spin), sigmaz(spin)
builder = OperatorBuilder(sys, auto_hermitian = true)
for site in l
    site_x = site + Bravais[1, 0]
    site_y = site + Bravais[0, 1]
    builder[site, site] = (1 + 0.1 * randn()) * sz
    builder[site, site_x] = (sz - im * sx) / 2
    builder[site, site_y] = (sz - im * sy) / 2

    builder[site, site] += 0.5 * E * site.x
    for site2 in adjacentsites(NearestNeighbor(l, 2), site)
        # note that the coefficient is 0.05, not 0.1, because every bond is counted twice
        builder[site, site2] += 0.05    
    end
end

c_site = l[!, x = 2, y = 2]
builder[c_site, c_site] += 0.5  # Added local potential on site (2, 2)

H2 = Hamiltonian(builder)
P = densitymatrix(H, info=false)
heatmap(localdensity(P))
```

Note that `builder[site1, site2]` is not a regular indexing operation. The return value of it is not an matrix or operator, but rather a special object that makes increment/decrement operations possible. Do not use this value in other contexts.

!!! tip
    Use `FastOperatorBuilder` instead of `OperatorBuilder` if you need to create a large Hamiltonian. It is a little 
    bit faster, but doesn't support direct assignment `builder[site, site] = ...`; use `+=` or `-=` instead.

## Gauge fields

There is a convenient way to add gauge fields to the Hamiltonian, and it is done by using the `GaugeField` objects.
It is an interface for different types of gauge fields, like the Landau gauge. To add field to the Hamiltonian, use the `field` keyword argument:

```@example 4
using LatticeModels, Plots
l = SquareLattice(4, 4)
H = tightbinding_hamiltonian(l, field=LandauGauge(0.1))
P = densitymatrix(H, info=false)
heatmap(localdensity(P))
```

This adds a magnetic field in the Landau gauge to the Hamiltonian: ``\overrightarrow{\mathcal{A}} = B x \overrightarrow{e_y}``. It adds a phase factor to the hopping terms, which is calculated using [Peierls substitution](https://en.wikipedia.org/wiki/Peierls_substitution).

Here are all types of gauge fields supported by this package:
- [`LandauGauge(B)`](@ref LandauGauge) — the Landau gauge, with the magnetic field `B` in the $z$ direction.
- [`SymmetricGauge(B)`](@ref SymmetricGauge) — the symmetric gauge, with the magnetic field `B` in the $z$ direction.
- [`PointFlux(Phi, center=(0, 0); gauge=:axial)`](@ref PointFlux) — a point flux in the point `center` with the flux `Phi`. The `gauge` argument can be either `:axial` (default) or `:singular`.
- [`GaugeField(f; n)`](@ref GaugeField) — a general magnetic field. The `f` is a function that takes a coordinate vector and returns the vector potential ``\mathcal{A}`` at this point. The line integrals are calculated using the `n`-point trapezoidal rule. Note that the `n` must be set explicitly.
- [`LineIntegralGaugeField(f)`](@ref LineIntegralGaugeField) — a general magnetic field. The `f` is a function that takes two coordinate vectors and returns the ``\int_{\vec{r}_1}^{\vec{r}_2} \mathcal{A} \cdot d\vec{r}`` line integral of the vector potential between these points.

You can pass these objects using the `field` keyword argument to the `tightbinding_hamiltonian`, `construct_hamiltonian`, and `OperatorBuilder` functions to add the gauge field to the Hamiltonian:

```@example 4
builder = OperatorBuilder(l, auto_hermitian = true, field = LandauGauge(0.1))
for site in l
    site_x = site + Bravais[1, 0]
    site_y = site + Bravais[0, 1]
    builder[site, site_x] = 1
    builder[site, site_y] = 1
end
H2 = Hamiltonian(builder)
H == H2
```

To find out more about operators, diagonalization and observables, proceed to the next section.