## Hamiltonian macro

A hamiltonian operator usually consists of a [diagonal part](@ref Diagonal-operators) and a [hopping part](@ref Hopping-operators). 
Sometimes a magnetic field is applied, in which case additional phase factors emerge according to [Peierls substitution](https://en.wikipedia.org/wiki/Peierls_substitution). 

Taking all this into account, the operator will now look something like this:

$$\hat{H} = 
\sum_i \hat{s}_i \hat{c}^\dagger_i \hat{c}_i + \left( \sum_{\text{adjacent }i, j} \hat{t}_{ij} \hat{c}^\dagger_j \hat{c}_i
\cdot e^{\frac{2\pi i}{\phi_0} \int_{r_i}^{r_j} \overrightarrow{A} \cdot \overrightarrow{dl}} + h. c. \right)$$

Such an operator can be easily constructed using the [`@hamiltonian`](@ref) macro. 
All you have to do is assign `lattice` and also `field` if needed, and then define the diagonal and hopping members using `@diag` and `@hop`.

Let's take a look at the example on the [Usage examples](@ref) page:

A simple Hubbard model hamiltonian for a square lattice is defined by the formula 
$$\hat{H} = \sum_\text{x-bonds} c^\dagger_i c_j + \sum_\text{y-bonds} c^\dagger_i c_j + h. c.$$

We can create a matrix for this operator on a `xsize`×`ysize` square lattice with the following code:

```@setup env
using LatticeModels
```

```@example env
Hubbard(xsize, ysize, field=NoField()) = @hamiltonian begin
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
\sum_i m_i c^\dagger_i \sigma_z c_i + \left(
\sum_\text{x-bonds} c^\dagger_i \frac{\sigma_z - i \sigma_x}{2} c_j + 
\sum_\text{y-bonds} c^\dagger_i \frac{\sigma_z - i \sigma_y}{2} c_j + 
h. c. \right)$$

Here we want the $m_i$ values to be represented by any number-typed `LatticeValue` - so the most convenient way to do it is to use tensor product notation:

```@example env
# The Pauli matrices
const σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]

# Initial hamiltonian: m=1 everywhere
ChernInsulator(m::LatticeValue{<:Number, :square}, field=NoField()) = @hamiltonian begin   
    lattice := lattice(m)
    field := field
    @diag m ⊗ σ[3]
    @hop (σ[3] - im * σ[1]) / 2 axis = 1
    @hop (σ[3] - im * σ[2]) / 2 axis = 2
end
nothing # hide
```

Note that the `field` parameter is defaulted to `NoField()` - this is a magnetic field object representing zero field. Other available magnetic field objects are:

- `LandauField(B::Real)` - Landau-calibrated uniform magnetic field
- `SymmetricField(B::Real)` - symmetrically calibrated uniform magnetic field
- `FluxField(Φ::Real, point::NTuple{2, Number})` - point magnetic field flux through `point`.

Most of these fields are designed for 2D lattices, but may also be applied to lattices with other dimensionality.

## Operator spectrum utilities

A `Spectrum` object contains eigenvalues and eigenvectors of some hermitian operator. 

It is a very convenient way to work with eigenvectors. Check this out:

```@repl env
sp = spectrum(Hubbard(5, 5))
sp[1]       # Get the first eigenstate (e. g. ground state)
sp[3]       # Get the third eigenstate
sp[E = 0]   # Get eigenstate with energy nearest to 0
sp[1:3]     # Create a new Spectrum with states from #1 to #3
E = eigvals(sp) # Get eigenvalues
sp_filled = sp[E .< 0]  # Create a new Spectrum with filled bands only
```

if you need an operator that projects onto some states, use [`projector`](@ref) or [`filled_projector`](@ref) functions:

```@repl env
P = projector(E -> E < 0, sp)
P = projector(sp_filled)    # The same thing
P = filled_projector(sp, 0) # The same thing
P = filled_projector(sp)    # The same thing
```