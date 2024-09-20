using Test, Documenter, LinearAlgebra, StaticArrays, Plots
using LatticeModels

include("test_workflows.jl")
include("test_lattice.jl")
include("test_field.jl")
include("test_operators.jl")
include("test_timedeps.jl")
include("test_currents.jl")

if VERSION >= v"1.10"
    doctest(LatticeModels)
end
