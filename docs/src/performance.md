This page is dedicated to various performance tips and tricks for `LatticeModels.jl`. Some of them are specific to the package, while others are general Julia tips.

## Add a site lookup table to the lattice

If a lattice is large, the site lookup can be slow. To speed up the process, one can add a lookup table to the lattice – this is a helper structure that maps the site coordinates to the site indices. This can be done the following way:

```julia
using LatticeModels
l = GrapheneRibbon(3000, 3000, boundaries=(:horizontal => true, :vertical => false))
ll = addlookuptable(l)
@time tightbinding_hamiltonian(l) # 4.6 s
@time tightbinding_hamiltonian(ll) # 3.0 seconds
```

As you can see, it starts making a difference for REALLY large lattices. Still, bear in mind that for a ``N``-site lattice the lookup table is constructed for ``\mathcal{O}(N)`` time, and the lookup of one site is ``\mathcal{O}(1)`` with the table and ``\mathcal{O}(\log N)`` without it, so eventually it pays off.

Note, however, that the result of `addlookuptable` is a new lattice. Moreover, this function is not type-stable, since the type inference is too complex for the compiler. Use it only outside of functions.

## Boost the `OperatorBuilder`

The `OperatorBuilder` is a very powerful tool for constructing operators on the lattice. However, when used the `default` way it is slower than `construct_operator`, especially for large lattices. The upside, however, is flexibility: the performance and memory consumption will be plausible in most cases for not very big lattices.

```julia
using LatticeModels
l = GrapheneRibbon(1000, 1000, boundaries=(:horizontal => true, :vertical => false))

function H1(l::AbstractLattice)
    builder = OperatorBuilder(l, auto_hermitian=true)
    for site in l
        builder[site, site] = 1.0
        builder[site, site + Bravais[1, 0]] = 0.2 - 0.1im
        builder[site, site + Bravais[0, -1]] = 0.4
        builder[site, site + Bravais[-1, 1]] = 0.6
    end
    return Hamiltonian(builder)
end
@time H1(l) # 1.82 s
```

If the Hamiltonian generation is a bottleneck (for example, for a large lattice or time-dependent evolution), you can adjust a more optimized way to construct the Hamiltonian:

```julia
function H2(l::AbstractLattice)
    builder = OperatorBuilder(UniformMatrixBuilder{ComplexF64}, l, auto_hermitian=true, col_hint=7)
    for site in l
        builder[site, site] = 1.0
        builder[site, site + Bravais[1, 0]] = 0.2 - 0.1im
        builder[site, site + Bravais[0, -1]] = 0.4
        builder[site, site + Bravais[-1, 1]] = 0.6
    end
    return Hamiltonian(builder)
end
@time H2(l) # 1.36 s
ll = addlookuptable(l) 
@time H2(ll) # 0.83 s - even faster with a lookup table
```

What happened here? The `UniformMatrixBuilder` is a more optimized way to construct a sparse matrix. It works well under the assumption that the number of nonzero elements in each column is approximately the same, hence the name. The `col_hint` parameter is a hint for the builder to preallocate the memory: we expect 7 nonzero elements in each column, one for the diagonal part and two per hopping.

Note that this exact algorithm is used under the hood of `construct_operator`:

```julia
@time construct_operator(ll, 1.0, 0.2 - 0.1im => Bravais[1, 0], 0.4 => Bravais[0, -1], 0.6 => Bravais[-1, 1]) # 0.43 s
```

As you can see, the `construct_operator` is still faster. Why? Because in our current version of `H2` the `site` site is looked up 5 times per iteration. However, in fact we already know its index, because it is a variable we iterate over. This can be fixed by using an integer index instead:

```julia
function H3(l::AbstractLattice)
    builder = OperatorBuilder(UniformMatrixBuilder{ComplexF64}, l, auto_hermitian=true, col_hint=7)
    i = 1
    for site in l
        builder[i, i] = 1.0
        builder[i, site + Bravais[1, 0]] = 0.2 - 0.1im
        builder[i, site + Bravais[0, -1]] = 0.4
        builder[i, site + Bravais[-1, 1]] = 0.6
        i += 1
    end
    return Hamiltonian(builder)
end
@time H3(ll) # 0.5 s - almost as good
```

By carefully implementing various techniques, we achieved a more than 3x speedup!

## Follow the Julia performance tips

You can find a lot of useful tips in the [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) section of the official documentation. Some of them are not directly related to this package, but can still be useful.

Most importantly:
- Write your code in functions and use them as much as possible. This allows the compiler to optimize the code.
- Avoid global variables, use function arguments instead. If global variables are really necessary in your case, use type annotations.
- Do not use abstract types in performance-critical code. Remember, `Integer`, `Real`, `Complex` are abstract types! 
  The concrete types meant here are respectively `Int`, `Float64`, `ComplexF64`.
- Make sure your functions are type-stable, e.g the type of each variable is well-defined at compile time. 
  Use `@code_warntype` to track down possible type instabilities.