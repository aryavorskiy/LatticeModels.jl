var documenterSearchIndex = {"docs":
[{"location":"currents/","page":"Currents","title":"Currents","text":"TODO","category":"page"},{"location":"library/#LatticeModels.jl","page":"Library","title":"LatticeModels.jl","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]","category":"page"},{"location":"library/#Lattice-creation","page":"Library","title":"Lattice creation","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"lattice.jl\"]","category":"page"},{"location":"library/#LatticeModels.Bravais","page":"Library","title":"LatticeModels.Bravais","text":"Bravais{N, NB}\n\nN-dimensional infinite Bravais lattice with NB sites in basis.\n\n\n\nBravais(translation_vectors[, basis])\n\nConstructs a Bravais lattice with given translation vectors and locations of basis sites relative to some unit cell. The basis argument can be omitted, in which case the lattice basis will consist of one site located in the bottom-left corner of the unit cell.\n\ntranslation_vectors argument must be an AbstractMatrix{<:Real} of size N×N, while basis must also be an  abstract matrix of size N×NB.\n\n\n\n\n\n","category":"type"},{"location":"library/#LatticeModels.Lattice","page":"Library","title":"LatticeModels.Lattice","text":"Lattice{LatticeSym, N, NB}\n\nA finite subset of a Brvais{N, NB}. LatticeSym is a Symbol which represents the type of the lattice (e. g. :square, :honeycomb). This makes Lattice object behavior known at compile-time, which allows to introduce various optimizations or to define specific plot recipes.\n\n\n\nLattice(sym, sz, bvs[, mask])\n\nConstructs a finite Lattice{sym, N, NB} as a subset of the bvs Bravais lattice. sz is a NTuple{N, Int} which represents how many times the unit cell of bvs was translated by each axis - these sites form a macro cell. mask, if defined, is a Vector{Bool} storing information about which of the sites from the macro cell are actually included in the lattice, and which are not.\n\nFor example, a 3×3 square lattice with its center site excluded is represented as Lattice(:square, (3, 3), Bravais([1 0; 0 1]), Bool[1, 1, 1, 1, 0, 1, 1, 1, 1])\n\nTo define a new type of lattice, create an alias for Lattice{YourSym, YourN, YourNB}. Refer to the docs for detailed explanation.\n\n\n\n\n\n","category":"type"},{"location":"library/#LatticeModels.LatticeSite","page":"Library","title":"LatticeModels.LatticeSite","text":"LatticeSite{N}\n\nA site of a Lattice{LatticeSym, N, NB} lattice. Fields:\n\nunit_cell: a set of translations along all axes representing the unit cell the site is located in.\nbasis_index: the number of site in the lattice basis.\n\nThis type is used to iterate over all sites of a Lattice{LatticeSym, N, NB}. The exact location of a LatticeSite can be found using the coords(lattice, site) function.\n\n\n\n\n\n","category":"type"},{"location":"library/#LatticeModels.SquareLattice","page":"Library","title":"LatticeModels.SquareLattice","text":"SquareLattice{N}\n\nBasically the same as Lattice{:square,N,1}.\n\n\n\nSquareLattice(sz::Int...)\n\nConstructs a square lattice with size sz.\n\n\n\n\n\n","category":"type"},{"location":"library/#LatticeModels.coords-Tuple{Lattice, LatticeSite}","page":"Library","title":"LatticeModels.coords","text":"coords(lattice::Lattice, site::LatticeSite) -> vector\n\nFinds the location in space of lattice site site on lattice lattice.\n\n\n\n\n\n","category":"method"},{"location":"library/#LatticeModels.radius_vector-Tuple{Lattice, LatticeSite, LatticeSite}","page":"Library","title":"LatticeModels.radius_vector","text":"radius_vector(lattice::Lattice, site1::LatticeSite, site2::LatticeSite) -> vector\n\nFinds the vector between two sites on a lattice according to possibly periodic boundary conditions (site2 will be translated along the macro cell to minimize the distance between them).\n\n\n\n\n\n","category":"method"},{"location":"library/#LatticeModels.sublattice-Union{Tuple{LatticeSym}, Tuple{Function, Lattice{LatticeSym}}} where LatticeSym","page":"Library","title":"LatticeModels.sublattice","text":"sublattice(lf::Function, l::Lattice) -> Lattice\n\nGenerates a a subset of lattice l by applying the lf function to its sites. The lf function must accept two positional arguments (a LatticeSite and a vector with its coordinates) and return a boolean value.\n\n\n\n\n\n","category":"method"},{"location":"library/#Lattice-values","page":"Library","title":"Lattice values","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"lattice_value.jl\"]","category":"page"},{"location":"library/#Lattice-operators","page":"Library","title":"Lattice operators","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"lattice_operator.jl\"]","category":"page"},{"location":"library/#Hoppings","page":"Library","title":"Hoppings","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"hoppings.jl\"]","category":"page"},{"location":"library/#Magnetic-fields","page":"Library","title":"Magnetic fields","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"field.jl\"]","category":"page"},{"location":"library/#Hamiltonian-tools","page":"Library","title":"Hamiltonian tools","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"hamiltonian.jl\"]","category":"page"},{"location":"library/#Currents","page":"Library","title":"Currents","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"currents.jl\"]","category":"page"},{"location":"library/#Unitary-evolution","page":"Library","title":"Unitary evolution","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules = [LatticeModels]\nPages   = [\"evolution.jl\"]","category":"page"},{"location":"library/#LatticeModels.evolution_operator-Tuple{Any, Real}","page":"Library","title":"LatticeModels.evolution_operator","text":"evolution_operator(H, t[, k])\n\nCalculates the unitary evolution operator using the formula\n\n$ \\mathcal{U}(t) = e^{-\\frac{1}{i\\hbar} \\hat{H} t} $\n\nArguments\n\nH: the hamiltonian matrix\nt: the evolution time\nk: if provided, the exponent will be calculated using a Taylor series expansion with order k\n\n\n\n\n\n","category":"method"},{"location":"library/#LatticeModels.@evolution-Tuple","page":"Library","title":"LatticeModels.@evolution","text":"@evolution [kwargs...] {rules...} for_loop\n\nGenerates an environment with defined hamiltonian and density matrices that evolve by certain laws. See Unitary evolution for more details.\n\n\n\n\n\n","category":"macro"},{"location":"operator/","page":"Lattice values and operators","title":"Lattice values and operators","text":"TODO","category":"page"},{"location":"plot/","page":"Visualization","title":"Visualization","text":"TODO","category":"page"},{"location":"advanced/","page":"Advanced options","title":"Advanced options","text":"TODO","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"TODO","category":"page"},{"location":"lattice/","page":"Defining a lattice","title":"Defining a lattice","text":"TODO","category":"page"},{"location":"#LatticeModels.jl","page":"Home","title":"LatticeModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides a set of tools to simulate different quantum lattice systems.","category":"page"},{"location":"#Package-features","page":"Home","title":"Package features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Bravais lattices with arbitrary geometry and any possible count of internal states on one sites.\nVersatile hamiltonian generation tools.\nBackend-independent computations: linear operators can be of any array type, allowing to use sparse or GPU arrays when needed.\nSmart unitary evolution macro reducing excessive computations where possible.\nPlots.jl integration.","category":"page"},{"location":"#Example-usage","page":"Home","title":"Example usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using LatticeModels\nusing LinearAlgebra, Plots\n\nl = SquareLattice(10, 10)\n\n# Define a Hofstadter model hamiltonian\nh(B) = @hamiltonian begin   \n    lattice := l\n    # Add hoppings along axis x and y\n    @hop axis = 1\n    @hop axis = 2\n    # Add magnetic field through (0, 0) point\n    field := FluxField(B, (0, 0))\nend\n\n# Calculate eigenvalues and eigenvectors\nsp = spectrum(h(0))\n\n# Find density matrix for filled bands (e. g. energy < 0)\nP_0 = filled_projector(sp)\n# Perform unitary evolution\nτ = 10\na = Animation()\n@evolution {\n    H := h(0.1 * min(t, τ) / τ)\n    P_0 --> H --> P\n} for t in 0:0.1:2τ\n    # Find the partial trace and plot it\n    density = diag_aggregate(m -> real(tr(m)), P)\n    plot(density, clims=(0,1))\n\n    # Show currents on the plot\n    plot!(DensityCurrents(H, P), arrows_scale=7)\n\n    title!(\"t = $t\")\n    frame(a)\nend\n\ngif(a, \"animation.gif\")","category":"page"}]
}
