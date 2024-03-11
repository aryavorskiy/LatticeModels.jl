This is a short (and not full) list of breaking changes introduced since v0.20:

- The `coordoperators` function now accepts a `Basis` instead of an integer as the second argument. You can use `Genericbasis(N)` instead of simply `N` if you do not want to specify the basis, but this is highly unrecommended.
- The `filled_projector` function is replaced with `densitymatrix`. The chemical potential `μ`, temperature `T` and particle statistics `statistics` can be set via keyword arguments. The default values are extracted from the hamiltonian, see docstring for more info.
- The `@hamiltonian` macro was removed. Use `construct_hamiltonian` or `OperatorBuilder`.
- The `@on_lattice` macro was removed, too, since the `LatticeArray`s were discontinued in favor of `Operator`s from QuantumOptics.jl. Common linear algebra operations are expected to work as usual.
- The `@field_def` macro was also removed. Subtype the `AbstractField` type explicitly or use `GaugeField` if you do not need performance.
- Many of the functions were renamed according to the style guide: `lattice_density` -> `localdensity`, `coord_values` -> `coordvalues` etc.
- The gauge field structs were renamed: `LandauField` -> `LandauGauge`, `SymmetricField` -> `SymmetricGauge`, `FluxField` -> `PointFlux`.
- The `@evolution` macro was removed. Use `Evolution` struct directly. 
- ...