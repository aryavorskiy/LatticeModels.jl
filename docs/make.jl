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
        "Tutorial" => "tutorial.md",
        "Step-by-step tutorial" => [
            "Defining the system" => "system.md",
            "Operators" => "operators.md",
            "Processing data" => "processing.md",
            "Currents" => "currents.md",
            "Unitary evolution" => "evolution.md"
        ],
        "Internals" => "internals.md",
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/aryavorskiy/LatticeModels.jl.git",
)
