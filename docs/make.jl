using Documenter
using LatticeModels

push!(LOAD_PATH,"../src/")

format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [joinpath("assets", "favicon.ico")]
)

makedocs(
    modules = [LatticeModels],
    sitename = "LatticeModels.jl",
    authors = "A. Yavorsky",
    format = format,
    pages = [
        "Home" => "index.md",
        "Tutorial" => [
            "Defining a lattice" => "lattice.md",
            "Lattice values" => "lattice_values.md",
            "Lattice operators" => "lattice_operator.md",
            "Hamiltonian generation" => "hamiltonian.md",
            "Currents" => "currents.md",
            "Unitary evolution" => "evolution.md"
        ],
        "Advanced options" => "advanced.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/aryavorskiy/LatticeModels.jl.git",
)
