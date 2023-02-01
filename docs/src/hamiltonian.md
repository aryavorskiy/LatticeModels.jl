## Hamiltonian macro

A hamiltonian operator usually consists of a [diagonal part](@ref Diagonal-operators) and a [hopping part](@ref Hopping-operators). 
Sometimes a magnetic field is applied, in which case additional phase factors emerge according to [Peierls substitution](https://en.wikipedia.org/wiki/Peierls_substitution). 

Taking all this into account, the operator will now look something like this:

$$\hat{H} = 
\sum_i^\text{sites} \hat{s}_i \hat{c}^\dagger_i \hat{c}_i + \left( \sum_{i, j}^\text{adjacent} \hat{t}_{ij} \hat{c}^\dagger_j \hat{c}_i
\cdot e^{\frac{2\pi i}{\phi_0} \int_{r_i}^{r_j} \overrightarrow{A} \cdot \overrightarrow{dl}} + h. c. \right)$$

Such an operator can be easily constructed using the [`@hamiltonian`](@ref) macro. 
All you have to do is assign `lattice` and also `field` if needed, and then define the diagonal and hopping members using `@diag` and `@hop`.

Let's take a look at the example on the [Usage examples](@ref) page:

A simple tight-binding model hamiltonian for a square lattice is defined by the formula 
$$\hat{H} = \sum_i^\text{sites} \left( c^\dagger_{i + \hat{x}} c_i + c^\dagger_{i + \hat{y}} c_i + h. c. \right)$$

We can create a matrix for this operator on a `xsize`×`ysize` square lattice with the following code:

```@setup env
using LatticeModels
```

```@example env
TightBinding(xsize, ysize, field=NoField()) = @hamiltonian begin
    lattice := SquareLattice(xsize, ysize)
    field := field
    @hop axis = 1   # x-bonds
    @hop axis = 2   # y-bonds
end
nothing # hide
```

Note that the keyword arguments for hopping operators are written as if they were passed to the [`hopping`](@ref) function.

For a Chern insulator the hamiltonian looks like this:

$$\hat{H} = 
\sum_i^\text{sites} m_i c^\dagger_i \sigma_z c_i + 
\sum_i^\text{sites} \left( 
c^\dagger_{i + \hat{x}} \frac{\sigma_z - i \sigma_x}{2} c_i + 
c^\dagger_{i + \hat{y}} \frac{\sigma_z - i \sigma_y}{2} c_i + 
h. c. \right)$$

Here we want the $m_i$ values to be represented by any number-typed `LatticeValue`, so the most convenient way to do it is to use tensor product notation:

```@example env
# The Pauli matrices
const σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]

# Initial hamiltonian: m=1 everywhere
ChernInsulator(m::LatticeValue{<:Number, :square}, field=NoField()) = @hamiltonian begin   
    lattice := lattice(m)
    field := field
    dims_internal := 2
    @diag m ⊗ σ[3]
    @hop (σ[3] - im * σ[1]) / 2 axis = 1
    @hop (σ[3] - im * σ[2]) / 2 axis = 2
end
nothing # hide
```

Here we must explicitly set the internal phase space dimension count via `dims_internal := 2`.

Some default hamiltonian formulas are already implemented - see [`TightBinding`](@ref), [`SpinTightBinding`](@ref), [`Haldane`](@ref) docstrings for more info.

!!! tip
    It is possible to set the matrix type of the hamiltonian operator at the generation time - use `arrtype := <Preferred type>` notation, this can be used e. g to reduce memory usage by switching to sparse arrays.
    Note however that some array types may not support this feature, in which case use the [`@on_lattice`](@ref) macro.

Note that the `field` parameter is defaulted to `NoField()` (this is a magnetic field object representing zero field). Other available magnetic field objects are:

- `LandauField(B::Real)` - Landau-calibrated uniform magnetic field
- `SymmetricField(B::Real)` - symmetrically calibrated uniform magnetic field
- `FluxField(Φ::Real, point::NTuple{2, Number})` - point magnetic field flux through `point`.

Most of these fields are designed for 2D lattices, but may also be applied to lattices with other dimensionality. You can define your own field types, see [Custom magnetic fields](@ref Custom-magnetic-fields).

## Operator spectrum utilities

A `Spectrum` object contains eigenvalues and eigenvectors of some hermitian operator. 

It is a very convenient way to work with eigenvectors. Check this out:

```@repl env
sp = spectrum(TightBinding(5, 5))
sp[1]       # Get the first eigenstate (e. g. the state with lowest energy)
sp[3]       # Get the third eigenstate
sp[E = 0]   # Get eigenstate with energy nearest to 0
sp[1:3]     # Create a new Spectrum with states from #1 to #3
E = eigvals(sp) # Get eigenvalues
sp_filled = sp[E .< 0]  # Create a new Spectrum with filled bands only
```

if you need an operator that projects onto some states, use [`projector`](@ref) or [`filled_projector`](@ref) functions:

```@repl env
P = projector(E -> E < 0, sp)
P = projector(sp_filled)                # The same thing
P = filled_projector(sp, 0)             # The same thing
P = filled_projector(sp)                # The same thing
```

Pass a lambda as a first argument to `filled_projector` to induce arbitrary state density. 
The lambda must take the energy and return the density as a real number.
The [`fermi_dirac`](@ref) and [`bose_einstein`](@ref) helper functions will generate such lambdas for the respective statistics.

```@repl env
PT = projector(E -> 1 / (exp(E) + 1), sp)
PT = projector(fermi_dirac(0, 1), sp)   # The same thing
```

!!! warning
    The [`spectrum`](@ref) function finds the eigenvalues and eigenvectors using `LinearAlgebra.eigen`.
    You will probably have to define it for array types that do not support default `LinearAlgebra` routines.

One more thing you can do with `Spectrum`s is finding the density of states (DOS) and the local density of states (LDOS). 
The DOS is the imaginary part of $\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ operator, whereas the LDOS is its partial trace.

It's highly recommended to use the built-in [`dos`](@ref) and [`ldos`](@ref) functions, because the `Spectrum` object stores the $\hat{H}$ operator already diagonalized, which makes calculations much faster.

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
