## Custom magnetic fields

You can induce magnetic field by setting the vector potential $\overrightarrow{A}$ in every point of the space.
This can be done by creating an `AbstractField` object: this is a lazy object that stores information about the magnetic field
and implements the [`LatticeModels.vector_potential`](@ref) function that calculates the above-mentioned vector potential.

The $\int_{r_1}^{r_2} \overrightarrow{A} \cdot \overrightarrow{dl}$ integral is found automatically with the [`LatticeModels.path_integral`](@ref) function.
The number of steps can be adjusted, but it is recommended to redefine this method for each new magnetic field type using an exact formula: this will improve performance and accuracy in most cases.

You can use a convenience macro to simplify definition of new field types:

```julia
@field_def struct LandauField(B::Number)
    # Define the vector potential function
    vector_potential(x) = (0, x*B)

    # Redefine the integrating function
    path_integral(p1, p2) = ((p1[1] + p2[1]) / 2) * (p2[2] - p1[2]) * B
end
```

Let's see what happened here. 
This macro created a new struct `LandauField <: AbstractField` with a constructor `LandauField(B::Number)`.

It also defined a suitable `LatticeModels.vector_potential` function. The only parameter in this function definition is `x`, which will be the X coordinate of the point (all other coordinate values will not be passed). 
The return type must be `Tuple` or `SVector` to ensure that it can be converted to `SVector` by the default `LatticeModels.path_integral` implementation.

This field object is compatible with lattices of any dimension count. Undefined elements of the $\overrightarrow{A}$ vector will be set to zero, and "excessive" ones will be discarded.

!!! tip
    To handle dimension-dependent behavior, you can use `vector_potential(p...) = (0, p[1]*B)` notation.
    Here `p` is a `SVector` which will allow you to yield different values depending on the dimension count.

You may notice that here the `LatticeModels.path_integral` function was also redefined. It simply takes two `SVector`s describing the coordinates of $p_1$ and $p_2$ and returns the value of the $\int_{p_1}^{p_2} \overrightarrow{A} \cdot \overrightarrow{dl}$ path integral between them.

!!! tip
    To improve the accuracy of the integration without redefining `LatticeModels.path_integral`, 
    set the default number of integration steps by adding `n_steps := <desired number>` to the struct definition.

## Custom array backends

You can freely change the type of underlying arrays in `LatticeArray`s by using the [`@on_lattice`](@ref) macro.
It will convert the internal array to the desired type and wrap it with the same `Basis`:

```@repl
using LatticeModels # hide
X, Y = coord_operators(SquareLattice(50, 50), 2)
XY = X * Y
using SparseArrays
sp_XY = @on_lattice sparse(XY)
```