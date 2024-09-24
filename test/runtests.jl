using LinearAlgebra, StaticArrays, Plots
using Test, Documenter
using LatticeModels

include("test_workflows.jl")
include("test_lattice.jl")
include("test_field.jl")
include("test_operators.jl")
include("test_timedeps.jl")
include("test_currents.jl")

if VERSION >= v"1.10"
    include("test_jet.jl")
    doctest(LatticeModels)
end
