## Lattice operators

```@setup env
using LatticeModels
```

A [`LatticeArray`](@ref LatticeModels.LatticeArray) is a wrapper type that maps an operator matrix or a wave function vector to a basis it is defined on.

A `Basis`, in turn, is generally a pair of a `Lattice` and an `Int` which is the number of internal states on every lattice site.

You can create `LatticeArray` representations of coordinate operators by using the `coord_operators` function. 

```@repl env
l = SquareLattice(5, 5)
bas = Basis(l, 2)
X, Y = coord_operators(bas) # Same as coord_operators(l, 2)
```

You can use arithmetic operators with `LatticeArray`s freely: after checking that the bases match they will be applied to the contained arrays,
and the result will be wrapped into a new `LatticeArray`. 
However, if you want to apply some other array operations (like `exp`), you can use the `@on_lattice` macro which will wrap the call as needed.

```jldoctest; setup=:(using LatticeModels, LinearAlgebra)
julia> X, Y = coord_operators(SquareLattice(5, 5), 2);

julia> X * sin(Y)
ERROR: MethodError: no method matching sin

julia> @on_lattice X * sin(Y)
50×50 LatticeOperator with inner type Matrix{Float64}
on Basis with 2-dimensional internal phase space
on 25-site square lattice on 5×5 macro cell
```

## Diagonal operators

A *lattice-diagonal* operator can be represented as 
$\hat{A} = \sum_i \hat{s}_i \hat{c}^\dagger_i \hat{c}_i$ , where $\hat{s}_i$ is an operator defined on the internal space and $\hat{c}^\dagger_i, \hat{c}_i$ are standard quantization operators.

Such an operator can be defined generally in two ways: 

**As a tensor product**

The $\hat{B}$ operator, being diagonal in the site basis, can easily be represented as a `LatticeValue`, while the $\hat{C}$ operator can be defined by a matrix.

```@repl env
l = HoneycombLattice(5, 5);
x, y = coord_values(l);
X, Y = coord_operators(l, 2);
op = (@. x + y) ⊗ [1 0; 0 1]
X + Y == op
```

**As a function call**

The `diag_operator` function is a low-level way to create a site-diagonal `LatticeOperator`. 
It accepts two arguments: one is a `Lattice` or a `Basis`, and the other is a function that can be written in the do-syntax same way as before - it must return a matrix which will affect the internal degrees of freedom on the site.

There are more possible ways to use this function, refer to [`diag_operator`](@ref) docstrings to find out more.

```@repl env
op1 = diag_operator(l) do site, (x, y)
    (x + y) * [1 0; 0 1]
end
op2 = diag_operator(Basis(l, 2), x .+ y)
X + Y == op1 == op2
```

Note that you can appply a function to every lattice-diagonal block of the operator matrix to create a `LatticeValue` with the `diag_reduce` function. 
For example, you can find a partial trace like this:

```@repl env
using LinearAlgebra
lv = diag_reduce(tr, op) # same as site_density(op)
lv == 2 .* x
```

[`site_density`](@ref) function is also applicable to wave functions represented as `LatticeVector`s.

## Hopping operators

A *hopping* operator is described by the following formula:

$$\hat{A} = \sum_{\text{adjacent }(i, j)} \hat{h} \hat{c}^\dagger_j \hat{c}_i + h. c.$$

The `Hopping` struct is a lazy representation of such an operator: it stores a matrix of the $\hat{h}$ operator 
and information about which sites should be considered *adjacent*. 
The relative location of the connected sites is defined by their respective indices in the lattice basis 
and the relative location of their unit cells as translations along lattice axes (represented as an integer vector).

A `Hopping` object can be created with the `hopping` function. Its only positional argument is 
the $\hat{h}$ operator matrix and keyword arguments `site_indices` and `translate_uc`. 

```@repl env
l = HoneycombLattice(5, 5)
hop1 = hopping(site_indices=(2, 1), translate_uc=[1, 0])
```

If the unit-cell translation is a single move along some lattice axis, you can instead set the `axis` keyword to the number of the translation axis. 
If there is no unit cell translation between the two sites, this argument can be omitted.

```@repl env
hop2 = hopping(site_indices=(2, 1), axis=2) # same as translate_uc=[0, 1]
hop3 = hopping(site_indices=(2, 1))
```

Note that no operator matrix argument was passed in the example above, which means that it was set to a 1×1 eye matrix automatically.

To generate a `LatticeOperator`, invoke the [`hopping_operator`](@ref) function with a `Lattice` and a `Hopping`:

```@repl env
hop_op = hopping_operator(l, hop1) + 
         hopping_operator(l, hop2) + 
         hopping_operator(l, hop3)
```

The bonds defined by the operator can be displayed on a plot in the following way:

```@example env
bs = bonds(hop_op)

using Plots
plot(bs)
plot!(l, show_excluded_sites=false, show_indices=false)
```

## Pair selectors

You can remove some of the hoppings by using a *pair selector*:
a *pair selector* is a lambda-like object which accepts a `Lattice` and two `LatticeSite`s to return whether the site pair is *selected* or not.
Passing such an object to the [`hopping_operator`](@ref) function as the first argument will remove the hoppings connecting *non-selected* pairs.

You can always write any lambda that accepts a `Lattice` and two `LatticeSite`s with similar aim to achieve similar results - but it's safer to use the [`PairSelector`](@ref LatticeModels.AbstractPairSelector) API, because this allows to define a selector for a greater lattice and use it for all sublattices safely.

For example, let's consider we want to split the lattice into two domains which are not connected to each other.
In such case we must firstly generate a `LatticeValue` defining which site corresponds to which domain:

```@example env
domains = @. abs(x) < 1 && abs(y) < 1 # A 2x2 square in the center of the lattice
selector = DomainsSelector(domains)
hop_op2 = hopping_operator(selector, l, hop1) + 
          hopping_operator(selector, l, hop2) + 
          hopping_operator(selector, l, hop3)

plot(bonds(hop_op2))
plot!(domains, cbar=false)
```

See also [`DomainsSelector`](@ref), [`PairLhsSelector`](@ref), [`PairLhsSelector`](@ref).