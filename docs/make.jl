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
    repo = Remotes.GitHub("aryavorskiy", "LatticeModels.jl"),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Manual" => [
            "Defining the lattice" => "manual/lattice.md",
            "Working with data" => "manual/latticevalue.md",
            "Constructing the Hamiltonian" => "manual/hamiltonian.md",
            "Operators" => "manual/operators.md",
            "Green's function" => "manual/greenfunction.md",
            "Currents" => "manual/currents.md",
            "Unitary evolution" => "manual/evolution.md"
        ],
        "Internals" => "internals.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/aryavorskiy/LatticeModels.jl.git",
)
