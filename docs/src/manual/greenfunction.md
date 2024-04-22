# Green's function

This chapter is dedicated to the (time-ordered) Green's function formalism. It is a powerful tool for studying many-body systems, especially in the context of condensed matter physics. The Green's function is a matrix-valued function of two variables, which can be used to calculate various physical quantities, such as the density of states or some topological invariants.

## Introduction

In the scope of this package, the (time-ordered) Green's function is defined by this formula:

```math
G_{\alpha\beta}(\omega) = \langle 0 | \hat{a}_{\alpha} \frac{1}{\omega - (\hat{H} - E_0)} \hat{a}_{\beta}^{\dagger} | 0 \rangle + q 
\langle 0 | \hat{a}_{\alpha}^{\dagger} \frac{1}{\omega + (\hat{H} - E_0)} \hat{a}_{\beta} | 0 \rangle, 
\hspace{1cm} 
q = \begin{cases} 1 & \text{for fermions} \\ -1 & \text{for bosons} \end{cases}
```

where ``H`` is the Hamiltonian of the system, ``| 0 \rangle`` is the ground state of the system (which is the vacuum for non-interacting systems), ``E_0`` is its energy and ``\hat{a}_{\alpha}`` and ``\hat{a}_{\beta}`` are the annihilation operators of the basis states. The Green's function is a matrix-valued function of the frequency ``\omega``.

This is how we can calculate the Green's function in the package:

```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(3, 2);

julia> H = tightbinding_hamiltonian(l);

julia> G = greenfunction(diagonalize(H))
Green's function for 6-site SquareLattice in 2D space

julia> G(0.1)       # Evaluate the Green's function at the frequency 0.1
Evaluated Green's function for 6-site SquareLattice in 2D space
6×6 Matrix{ComplexF64}:
 -0.209531+0.0im  0.0321264+0.0im  …  -0.108521+0.0im    1.04223+0.0im
 0.0321264+0.0im  -0.209531+0.0im       1.04223+0.0im  -0.108521+0.0im
  -1.05308+0.0im   0.212743+0.0im      -1.05308+0.0im   0.212743+0.0im
  0.212743+0.0im   -1.05308+0.0im      0.212743+0.0im   -1.05308+0.0im
 -0.108521+0.0im    1.04223+0.0im     -0.209531+0.0im  0.0321264+0.0im
   1.04223+0.0im  -0.108521+0.0im  …  0.0321264+0.0im  -0.209531+0.0im

julia> site1 = l[!, x=1, y=1]; site2 = l[!, x=2, y=2];

julia> G[site1, site2]
Green's function element (6 â†, 0 â bands)

julia> G[site1, site2](0.1) == G(0.1)[site1, site2]
true
```

Note that the Green's function is evaluated by exactly diagonalizing the Hamiltonian. This is not the most efficient way to calculate it, but it is the easiest to implement. Use brackets `G[site1, site2]` to get the Green's function between two sites, and call it with a frequency to evaluate it. 

There are three main types that provide the Green's function functionality:
- The `GreenFunction`, which is the main type for Green's functions. It stores the diagonalized Hamiltonian in an efficient way to optimize the Green's function calculation. Can be indexed with sites (or `(site, internal_ind)` tuples) and called with a frequency.
- The `GreenFunctionElement`, which is an element of the Green's function matrix. It stores the bands that correspond to the Green's function element. Can be called with a frequency.
- The `GreenFunctionEval`, which is the result of the Green's function evaluation. It stores the evaluated Green's function and the lattice it corresponds to. Can be indexed with sites.

Thus, `G[site1, site2]` creates a `GreenFunctionElement` object, while `G(0.1)` creates a `GreenFunctionEval` object. This explains why `G[site1, site2](0.1)` is the same as `G(0.1)[site1, site2]`, but the former is more efficient — it doesn't have to evaluate the whole Green's function matrix.

## Many-body Green's function

The Green's function can be evaluated for many-body systems as well. The formula remains the same, but the process of calculating it is a bit more complicated. You have to provide three Hamiltonian matrices for the `N-1`, `N` and `N+1` particle systems. Here is how it's done:

```@example 2
using LatticeModels, Plots
l = SquareLattice(4, 4)
H = bosehubbard(l, 2, U = 10) # Bose-Hubbard model, 2 particles
Hp1 = bosehubbard(l, 3, U = 10) # Bose-Hubbard model, 3 particles
Hm1 = bosehubbard(l, 1, U = 10) # Bose-Hubbard model, 1 particle
G = greenfunction(H, Hp1, Hm1, showprogress=false)
nothing # hide
```

Here is what happens in the code above:
1. Firstly, we find the ground state of the `N`-particle system and its energy. These will be ``| 0 \rangle`` and ``E_0`` in the Green's function formula.
2. Act with the annihilation and creation operators on ``| 0 \rangle``.
3. Use [`greenfunction`](@ref) to calculate the Green's function for the `N`-particle system.
   - The Green's function will be evaluated by exactly diagonalizing the Hamiltonian - this is why `greenfunction` requires the Hamiltonians for `N-1` and `N+1` particles.
   - For large Hamiltonian matrices, full diagonalization can be slow. In this case, the Lanczos algorithm with `` \hat{a}_{\alpha}^{\dagger} | 0 \rangle`` as initial vectors will be used. 

This is how you do it in code, line-by-line:

```@example 2
using LatticeModels, Plots
l = SquareLattice(4, 4)
H = bosehubbard(l, 2, U = 10) # Bose-Hubbard model, 2 particles
E, psi = findgroundstate(H)

Hm1 = bosehubbard(l, 1, U = 10) # Bose-Hubbard model, 1 particle
Hp1 = bosehubbard(l, 3, U = 10) # Bose-Hubbard model, 3 particles
G = greenfunction(psi, Hp1, Hm1, E0=E, showprogress=false)
nothing # hide
```

This is equivalent to the example before, but it allows you to set the ground state of the system manually.

Note that you can pass additional keyword arguments to `greenfunction` to control the diagonalization process. For example, the `routine` keyword argument allows you to choose the diagonalization routine, and the other keyword arguments are passed to the [`diagonalize`](@ref) function, see its documentation for more information.

!!! warning 
    Remember the order of the Hamiltonians in the `greenfunction` function. The first Hamiltonian should correspond to the `N`-particle system, the second to the `N+1`-particle system and the third to the `N-1`-particle system. If you pass them in the wrong order, an error will be thrown.

## Density of states

The density of states (DOS) is a physical quantity that can be calculated from the Green's function. It is defined as:

```math
DOS(\omega) = \frac{1}{\pi} \text{Im} \text{Tr} G(\omega - i\delta)
```

where ``\delta`` is the broadening. The DOS is a scalar function of the frequency ``\omega``. To showcase, let us calculate the DOS for the Bose-Hubbard model using the Green's function we calculated before:

```@example 2
df = dos(G, broaden=0.1)
plot(df, lab="", xlab="ω", ylab="DOS")
```

This will plot the DOS for the Bose-Hubbard model with a broadening of `0.1`. The DOS is a scalar function of the frequency, so it can be plotted as a line plot. 

Here `dos` produced a function that takes the frequency as an argument and returns the DOS at that frequency. You can also calculate the DOS for a specific frequency:

```@repl 2
dos(G, 0.1, broaden=0.1)
dos(G, 0.1, broaden=0.1) == df(0.1)
```

You can also calculate the local density of states (LDOS) using the [`ldos`](@ref) function. Here is how you do it:

```@example 2
ld = ldos(G, 0.1, broaden=0.1)
plot(ld, st=:shape)
```

Note that `ld` here is a [`LatticeValue`](@ref) object, which can be plotted as a shape plot. 

!!! tip
    To efficiently calculate the LDOS on one site, use this notation:

    ```@repl 2
    site = l[!, x=1, y=1]
    ld_value = dos(G[site, site], broaden=0.1)
    ld[site] == ld_value
    ```