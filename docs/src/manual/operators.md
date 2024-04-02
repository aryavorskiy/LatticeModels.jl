# States and Operators

In the previous chapter we have seen how to define a Hamiltonian for a lattice model. In this chapter we will see how to work with other types of operators, such as measurements and diagonalizing.

## Builtins

There are some basic building blocks that can be used to create states and operators. These are the `basisstate`, `diagonaloperator` and `transition` functions, taken from [QuantumOptics.jl](https://docs.qojulia.org/quantumobjects/operators/) package, but extended to work with lattices.

The `basisstate` function creates a state vector with a in the basis of a lattice. Normally it takes a `QuantumOptics.Basis` and an index as arguments, and returns a `Ket` object. However, you can use a lattice and a site as arguments to create a state vector in the basis of the lattice:

```@repl 0
using LatticeModels
l = SquareLattice(6, 6);
site = l[!, x = 3, y = 3]
l_bas, ind = LatticeBasis(l), site_index(l, site);
psi1 = basisstate(l_bas, ind);  # The pure QuantumOptics.jl way
psi2 = basisstate(l, site);     # The LatticeModels.jl way
psi1 == psi2
```

If you work in a composite system (e. g. a lattice with a spin degree of freedom), you can construct the state as tensor product of the states for each subsystem:

```@repl 0
spin = SpinBasis(1//2)
psi_composite = basisstate(spin, 1) ⊗ basisstate(l, site)   # Note the order!
```

!!! note
    The order of the tensor product is important. **The first argument is the on-site degrees of freedom, and the second is the lattice**. This convention is consistent in the rest of the package — the reason behind this is performance of [`construct_operator`](@ref) and [`OperatorBuilder`](@ref).

    If this order is not followed, you will probably get an error somewhere in your calculations.

The `diagonaloperator` function creates a diagonal operator in the basis of a lattice. Normally it takes a `QuantumOptics.Basis` and a vector of values (or a single value) as arguments, and returns an `Operator` object. There are several convenient ways to use this function with lattices. As an example let's consider the position operator in the basis of a lattice:

```@repl 0
xval = coordvalue(l, :x)
X_1 = diagonaloperator(l_bas, xval.values); # The pure QuantumOptics.jl way
X_2 = diagonaloperator(l, :x);              # The LatticeModels.jl way, with a site property
X_3 = diagonaloperator(xval);               # The LatticeModels.jl way, with a LatticeValue
X_1 == X_2 == X_3
```

This notation allows converting any `LatticeValue` or [Site parameter](@ref Sites) to an operator. Hence, `diagonaloperator(l, Coord(1))` is also valid and will return the same operator.

!!! tip
    Generally, a good value to create a custom diagonal operator is by using the `LatticeValue` approach. 
    
    Consider this example: in the Haldane model the diagonal part is ``m`` on the A sublattice and ``-m`` on the B sublattice. You can create this operator with the following code:

    ```julia
    l = HoneycombLattice(6, 6)
    m = 3
    ms = LatticeValue(l) do site
        site.index == 1 ? m : -m
    end
    Op = diagonaloperator(ms)
    ```

Also note the [`coordoperator`](@ref) and [`coordoperators`](@ref) functions that do the same thing as [`coordvalue`](@ref) and [`coordvalues`](@ref), but return the operator instead of the value:

```@repl 0
X, Y = coordoperators(l);
X == X_1
```

## Measurements

The most common type of measurements is the local density: the average number of particles at each site. This can be calculated using the [`localdensity`](@ref) function — it takes a state (a `QuantumOptics.Ket` vector or a `QuantumOptics.Operator` representing the density matrix) and returns a `LatticeValue` with the density at each site.

```@example 1
using LatticeModels, Plots
l = GrapheneRibbon(6, 4)
H = tightbinding_hamiltonian(l)
d = localdensity(groundstate(H))
plot(d)
```

The `localdensity` function uses the following formula to calculate the density at each site: ``\rho_i = \text{Tr}(\hat{n}_i \hat{\rho})``, where ``\hat{n}_i`` is the number operator at site ``i`` and ``\hat{\rho}`` is the density matrix. Note that if the values are complex, the function will return the real part of them. This is what makes the next example work.

The [local Chern marker](https://arxiv.org/abs/2304.13708) is a quantity that can be used to detect topological phases in a lattice model. It can be calculated using the following formula:

```math
\mathcal{C}(r) = 4\pi \text{Im} \langle r | P X P Y P | r \rangle
```

Here ``P`` is the projector onto the occupied states (e. g. the density matrix), and ``X`` and ``Y`` are the position operators. The local Chern marker is a real number that can be calculated for each site in the lattice.

Let's do this for a QWZ model Hamiltonian on a square lattice:

```@example 1
l = SquareLattice(6, 6)
sys = l ⊗ SpinBasis(1//2)
ms = ones(l)
ms[x = 3 .. 4, y = 3 .. 4] .= -1
H = qwz(l, ms)

P = densitymatrix(H, mu=0, statistics=FermiDirac)
X, Y = coordoperators(sys)
c = localdensity(-4π * im * P * X * P * Y * P)
heatmap(c, title="Local Chern marker")
```

## Diagonalizing

To diagonalize a Hamiltonian or any other operator, you can use the [`diagonalize`](@ref) function. It takes an operator and returns a `EigenSystem` object with the eigenvalues and eigenvectors of the operator.

```@example 2
using LatticeModels, Plots
l = GrapheneRibbon(6, 4)
H = haldane(l, 0.1, 1)
eig = diagonalize(H)
```

This struct simplifies the access to the eigenvalues and eigenvectors of the operator. You can access the eigenvalues with `eig.values`, and eigenvectors as `Ket`s can be obtained with the bracket notation `eig[i]` or `eig[value = E]`:

```@example 2
# The states are sorted by real part of the eigenvalues, so
psi = eig[1]                        # `psi` is the ground state
psi2 = eig[value = 0]               # `psi2` is the state with zero energy
p = plot(layout = @layout[a b; c], size=(800, 800))
plot!(p[1], localdensity(psi), title="Ground state")
plot!(p[2], localdensity(psi2), title="Zero energy state")
scatter!(p[3], eig.values, title="Spectrum", lab="")
```

!!! tip
    You can find the ground state of a Hamiltonian in one line using the [`groundstate`](@ref) function:

    ```julia
    psi = groundstate(H)
    ```

    To evaluate both the ground state and its energy, use the [`findgroundstate`](@ref) function:

    ```julia
    E0, psi = findgroundstate(H)
    ```

The `diagonalize` function under its hood uses the `eigen` function from the `LinearAlgebra` standard library. 
However, this does not work for non-trivial matrix types (e. g. sparse matrices or GPU arrays). For such cases you 
can pass the second argument to the `diagonalize` function, which is a `Symbol` indicating the method to use for 
diagonalization. To use the `eigsolve` function from the [`KrylovKit.jl`](https://jutho.github.io/KrylovKit.jl/stable/) 
package, you can pass `:krylovkit` as the second argument. Since this solves the eigenvalue problem iteratively, you 
can also pass the keyword arguments: `n` for the number of eigenvalues to compute, `v0` for the initial guess, and 
the keyword arguments for the `eigsolve` function.

```@example 2
l = SquareLattice(100, 100)             # A really big lattice
H = tightbinding_hamiltonian(l)
eig = diagonalize(H, :krylovkit, n=9)   # Compute only 9 eigenvalues with smallest real part
p = plot(layout=9, leg=false, size=(1000, 900))
for i in 1:9
    plot!(p[i], localdensity(eig[i]), title="E = $(round(eig.values[i], digits=5))",
        ms=2, msw=0, msa=0)             # Plot with small markers with no outline
end
plot!()
savefig("gs_density.png")   # hide
nothing                     # hide
```
![](gs_density.png)

`EigenSystem` objects have a wide range of applications in this package. One of them is creating equilibrium states for a given Hamiltonian. This can be done using the [`densitymatrix`](@ref) function, which is described in the next section.

## Density matrix

After you diagonalize a Hamiltonian, you can calculate the [density matrix](https://en.wikipedia.org/wiki/Density_matrix) for the system. Use the [`densitymatrix`](@ref) function to do this:

```@example 3
using LatticeModels
l = SquareLattice(6, 6)
sys = System(l, SpinBasis(1//2), mu=0, statistics=FermiDirac, T=0.1)
H = tightbinding_hamiltonian(sys)
eig = diagonalize(H)
P1 = densitymatrix(eig)                 # Use the default parameters from the `System`
P2 = densitymatrix(eig, 
    statistics=BoseEinstein, T=0, mu=1) # Or you can override them
nothing # hide
```

Note that the `densitymatrix` function can also be applied to a `Hamiltonian` object, in which case it will first diagonalize the Hamiltonian and then calculate the density matrix:

```@example 3
P1_1 = densitymatrix(H)
@assert P1 ≈ P1_1
P2_1 = densitymatrix(H, statistics=BoseEinstein, T=0, mu=1)
@assert P2 ≈ P2_1
nothing # hide
```

The `densitymatrix` function uses a simple formula to calculate the density matrix:

```math
\hat{\rho} = \sum_{i} \rho(E_i) | \psi_i \rangle \langle \psi_i |
```

Here ``E_i`` are the eigenvalues of the operator, ``\psi_i`` are the corresponding eigenvectors, and ``\rho(E)`` is the distribution function defined by the `statistics`, `T`, and `mu` parameters. By default, when no additional parameters are passed to the `System` or `densitymatrix`, the density matrix will be a thermal state at zero temperature.

The basis for these computations is the [`projector`](@ref) function, which takes a
function `p` and an `EigenSystem` object `d` that represents the diagonalized operator ``\hat{O}``. The return value is an operator ``\hat{P}`` defined by the formula:

```math
\hat{O} = \sum_{i} E_i | \psi_i \rangle \langle \psi_i |, \hspace{1cm}
\hat{P} = \sum_{i} p(E_i) | \psi_i \rangle \langle \psi_i |
```

Here ``E_i`` are the eigenvalues of the operator ``\hat{O}``, and ``| \psi_i \rangle`` are the corresponding eigenvectors. The function `p` is applied to the eigenvalues to obtain the diagonal elements of the density matrix. Here is an example of how to use this function:

```julia
l = SquareLattice(6, 6)
H = tightbinding_hamiltonian(l)
eig = diagonalize(H)
P1 = projector(x -> x < 0, eig)             # Projector onto the states with energy < 0
P2 = projector(x -> 1 / (1 + exp(x)), eig)  # Fermi-Dirac distribution
P3 = projector(eig[1:4])                    # Projector onto the first 4 states
# Note how we slice the `eig` object to get the first 4 states
```
