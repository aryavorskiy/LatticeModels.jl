using Documenter
using LatticeModels

push!(LOAD_PATH,"../src/")

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    assets = [joinpath("assets", "favicon.ico")],
    size_threshold_ignore = ["library.md"]
)

makedocs(
    modules = [LatticeModels],
    sitename = "LatticeModels.jl",
    authors = "A. Yavorsky",
    format = format,
    repo = Remotes.GitHub("aryavorskiy", "LatticeModels.jl"),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Manual" => [
            "Defining the lattice" => "manual/lattice.md",
            "Adding bonds" => "manual/bonds.md",
            "Working with 'LatticeValue's" => "manual/latticevalue.md",
            "Constructing the Hamiltonian" => "manual/hamiltonian.md",
            "States and Operators" => "manual/operators.md",
            "Green's function" => "manual/greenfunction.md",
            "Currents" => "manual/currents.md",
            "Evolution" => "manual/evolution.md"
        ],
        "API" => "library.md",
        "Internals" => "internals.md"
    ]
)

deploydocs(
    repo = "github.com/aryavorskiy/LatticeModels.jl.git",
)
