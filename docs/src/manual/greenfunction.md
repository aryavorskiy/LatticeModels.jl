# Green's function

This chapter is dedicated to the Green's function formalism. It is a powerful tool for studying many-body systems, especially in the context of condensed matter physics. The Green's function is a matrix-valued function of two variables, which can be used to calculate various physical quantities, such as the density of states or some topological invariants.

## Introduction

In the scope of this package, the Green's function is defined by this formula:

```math
G_{\alpha\beta}(\omega) = \langle 0 | \hat{a}_{\alpha} \frac{1}{\omega - (\hat{H} - E_0)} \hat{a}_{\beta}^{\dagger} | 0 \rangle - 
\langle 0 | \hat{a}_{\alpha}^{\dagger} \frac{1}{\omega + (\hat{H} - E_0)} \hat{a}_{\beta} | 0 \rangle
```

where ``H`` is the Hamiltonian of the system, ``| 0 \rangle`` is the ground state of the system (which is the vacuum for non-interacting systems), ``E_0`` is its energy and ``\hat{a}_{\alpha}`` and ``\hat{a}_{\beta}`` are the annihilation operators of the basis states. The Green's function is a matrix-valued function of the frequency ``\omega``. Note that this is the formula for bosons — fermions are not supported yet.

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

Note the three main types that provide the Green's function functionality:
- The `GreenFunction`, which is the main type for Green's functions. It stores the diagonalized Hamiltonian in an efficient way to optimize the Green's function calculation. Can be indexed with sites and called with a frequency.
- The `GreenFunctionElement`, which is an element of the Green's function matrix. It stores the bands that correspond to the Green's function element. Can be called with a frequency.
- The `GreenFunctionEval`, which is the result of the Green's function evaluation. It stores the evaluated Green's function and the lattice it corresponds to. Can be indexed with sites.

Thus, `G[site1, site2]` creates a `GreenFunctionElement` object, while `G(0.1)` creates a `GreenFunctionEval` object. This explains why `G[site1, site2](0.1)` is the same as `G(0.1)[site1, site2]`, but the former is more efficient — it doesn't have to evaluate the whole Green's function matrix.

## Many-body Green's function

## Density of states