In this chapter we will learn how to create Hamiltonians and other operators. There are several ways to do this in this package
 - from simple `tightbinding_hamiltonian`s to the flexible `build_operator` function to `OperatorBuilder`s.

## Simple tight-binding Hamiltonians

The simplest tight-binding Hamiltonian of possible is described by this formula:

$$\hat{H} = 
t \sum_{i, j}^\text{adjacent} \hat{b}^\dagger_j \hat{b}_i + 
t' \sum_{i, j}^\text{second adjacent} \hat{b}^\dagger_j \hat{b}_i + 
t'' \sum_{i, j}^\text{third adjacent} \hat{b}^\dagger_j \hat{b}_i + \ldots$$

Usually `t` is equal to `1`, while all other factors are zero.
Sometimes magnetic field is applied, in which case additional factors emerge according to [Peierls substitution](https://en.wikipedia.org/wiki/Peierls_substitution). 
This adds $\frac{2\pi}{\phi_0} \int_{r_i}^{r_j} \overrightarrow{A} \cdot \overrightarrow{dl}$ phase to every hopping.

To generate the Hamiltonian operator, use the [`tightbinding_hamiltonian`](@ref) function:

```julia
l = SquareLattice(10, 10, boundaries = ([10, 0] => true, [0, 10] => true))
H = tightbinding_hamiltonian(l)     # Yes, it is that simple
# t1 stands for t, t2 - for t', etc.
H2 = tightbinding_hamiltonian(l, t1 = 2, t2 = 1, t3 = 0)

s = System(l, N = 2, T = 1)
H3 = tightbinding_hamiltonian(s)    # You can use a `System` as well
```

Applying magnetic field is a bit more complicated. You can use a pre-defined magnetic field from the library or define your own. In the latter case 
you must provide the number of integration steps to calculate the phase in the Peierls substitution.

```julia
# 0.1 flux quantum per 1x1 plaquette, Landau gauge
H4 = tightbinding_hamiltonian(l, field=LandauField(0.1))

# Define the very same magnetic field
my_landau_field = MagneticField(n = 100) do crd
    x, y = crd
    # The function must return the vector potential as a `Tuple`
    return (0, 0.1 * x)
end
H5 = tightbinding_hamiltonian(l, field=my_landau_field)
```
The pre-defined magnetic field types include:
- `LandauField(B::Real)` - Landau gauge uniform magnetic field,
- `SymmetricField(B::Real)` - symmetric gauge uniform magnetic field,
- `FluxField(Φ::Real, point::NTuple{2, Number})` - point magnetic field flux through `point`.
  
Note that using a user-defined `MagneticField` will most probably result in slower operator generation. See [Custom magnetic fields](@ref) for more info.

## Other basic operators

To create spatial coordinate operators, use the [`coord_operators`](@ref) function. If you want one specific coordinate operator, use [`coord_operator`](@ref). These functions accept a `Sample` as first argument.

```julia
X = coord_operator(s, :x)       # X spatial coordinate operator
J1 = coord_operator(s, p"j1")   # Unit cell 1st coordinate operator
X, Y = coord_operators(s)       # Shorthand for generating both coordinate operators
X, Y = coord_operators(l)       # Also accepts a `Lattice` and optional internal basis
```

If more fine-grained control over the hoppings is required, you can make use of transition operators:

```julia
site1 = l[x = 2, y = 3]
site2 = l[x = 3, y = 2]
t12 = transition(s, site1, site2)
t34 = transition(s, 3, 4)           # 3 and 4 are integer indices of sites
H6 = H3 + t12 + t12' + t34 + t34'
```

## Builder function

The `build_hamiltonian` function allows to created a complicated hamiltonian in one function call. Let's get to know how it works.

Like `tightbinding_hamiltonian`, this function accepts a `System` (and all following default behavior) as first positional argument 
and magnetic field via `field` keyword. All other positional arguments describe different terms of the resulting hamiltonian. Each term is a `A => B` pair, where `A` describes an operator acting on the internal space, and `B` acts on the lattice space. The term is thus a tensor product of these operators.

Let's take a look at an example. We want to add spin-orbital interaction to our Hamiltonian. One way is to use
the QWZ model, whose Hamiltonian like this:

$$\hat{H} = 
\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i + 
\sum_i^\text{sites} \left( 
c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i + 
c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i + 
h. c. \right)$$

And this is how we can write a function that creates a hamiltonian in such a model.

```@setup env
using LatticeModels
```

```@example env
const σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]

qwzmodel(m::LatticeValue{<:Number, <:SquareLattice}; field=NoField()) = 
build_hamiltonian(lattice(m), SpinBasis(1//2), field=field,
    σ[3] => m,
    (σ[3] - im * σ[1]) / 2 => Bonds(axis = 1),
    (σ[3] - im * σ[2]) / 2 => Bonds(axis = 2))
nothing # hide
```

Let's explain what happens here line-by-line. 
Firstly, the parameter `m` is a `LatticeValue` - a struct which contains some value defined on the lattice's sites.
This is a very powerful tool, but all that we need now is that `LatticeValue{<:Number, <:SquareLattice}` stores numeric values on a square lattice.
We use implicit `Sample` construction here by passing the lattice `m` is defined on and a `SpinBasis(1//2)` (one half spin, obviously) as the internal 
dergee of freedom.

The `σ[3] => m` line describes the $\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i$ term. $m_i$ values are defined by the `LatticeValue` provided.

The next two lines describe the hopping terms in x- and y-direction. Take a closer look: `(σ[3] - im * σ[1]) / 2` on the left-hand side of the pair 
is the $\frac{\sigma_z - i \sigma_x}{2}$ operator acting on the spin of the particle. And `Bonds(axis = 1)` can be read as "hop from site to its neighbor along lattice axis #1". In context of a square lattice it is the same as the `x` coordinate axis.

We can use this function to find the ground state density matrix and plot the local Chern marker. Here we will also showcase 
part of the functionality `LatticeValue`s provide, see more in [Processing results](@ref).

```@example env
l = SquareLattice(10, 10)
m = ones(l)         # m_i will be equal to one everywhere...
m[x = 4..7] .= -1   # Except this ribbon
H = qwzmodel(m)

X, Y = coord_operators(l, SpinBasis(1//2))  # Coordinate operators
P = densitymatrix(H)                        # Density matrix
C = 4pi * im * P * X * P * Y * P
lcm = lattice_density(C)                    # Local Chern marker

using Plots
plot(lcm)
```

!!! note
    [`build_operator`](@ref) is actually the same function as [`build_hamiltonian`](@ref), but its output is an `Operator` rather than a `Hamiltonian`. 
    The only difference is that a `Hamiltonian` stores information about the `System` it is defined on. It is used to build a
    density matrix in a more convenient way. 

## The `OperatorBuilder`

!!! info "Docs under construction"
    This manual page is not written yet. I'm working on that. You can watch [a video on YouTube](https://www.youtube.com/watch?v=dQw4w9WgXcQ) meanwhile.

## Diagonalizing operators

!!! danger "Outdated docs"
    All information below is outdated and does not refer to the current functionality of the package.

A `Eigensystem` object contains eigenvalues and eigenvectors of some hermitian operator. 

It is a very convenient way to work with eigenvectors. Check this out:

```@repl env
dg = diagonalize(tightbinding_hamiltonian(SquareLattice(5, 5)))
dg[1]       # Get the first eigenstate (e. g. the state with lowest energy)
dg[3]       # Get the third eigenstate
dg[value = 0]   # Get eigenstate with energy nearest to 0
dg[1:3]     # Create a new Eigensystem with states from #1 to #3
E = dg.values   # Get eigenvalues
diag_filled = dg[E .< 0]    # Create a new Eigensystem with filled bands only
```

You can use [`densitymatrix`](@ref) or [`projector`](@ref) (only in zero-temperature cases) functions to find a density matrix of a system.

```julia
P = projector(E -> E < 0, dg)
P = projector(diag_filled)              # The same thing
P = densitymatrix(dg, μ = 0, T = 0)     # The same thing
P = densitymatrix(dg)                   # The same thing
```

!!! warning
    The [`diagonalize`](@ref) function finds the eigenvalues and eigenvectors using `LinearAlgebra.eigen`.
    You will probably have to define it for array types that do not support default `LinearAlgebra` routines.

One more thing you can do with `Eigensystem`s is finding the density of states (DOS) and the local density of states (LDOS). 
The DOS is the imaginary part of $\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ operator, whereas the LDOS is its partial trace.

It's highly recommended to use the built-in [`dos`](@ref) and [`ldos`](@ref) functions, because the `Eigensystem` object stores the $\hat{H}$ operator already diagonalized, which makes calculations much faster.

```@example env
using Plots

p = plot(layout=@layout[_ a{0.5w} _; grid(1, 2)], size=(800, 800))
plot!(p[1], -4:0.1:4, dos(sp, 0.2), title="DOS with δ = 0.2")
plot!(p[2], ldos(sp, 1, 0.2), title="LDOS at E = 1 with δ = 0.2")

# These two lines produce exactly the same result as the previous
ldos_fun = ldos(sp, 0.2)
plot!(p[3], ldos_fun(1), title="LDOS at E = 1 with δ = 0.2 (via function)")
```

!!! tip
    If you want to find the LDOS for the same `spectrum` and $\delta$, but for many different energy values, consider using `ldos(spectrum, delta)`: it improves performance dramatically.
