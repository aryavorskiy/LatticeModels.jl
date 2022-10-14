using Documenter
using LatticeModels

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [LatticeModels],
    sitename = "LatticeModels.jl",
    authors = "A. Yavorsky",
    pages = [
        "Home" => "index.md",
        "Tutorial" => [
            "Defining a lattice" => "lattice.md",
            "Lattice values" => "lattice_values.md",
            "Hamiltonian generation" => "hamiltonian.md",
            "Currents" => "currents.md",
            "Unitary evolution" => "evolution.md",
            "Visualization" => "plot.md"
        ],
        "Advanced options" => "advanced.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/aryavorskiy/LatticeModels.jl.git",
)
