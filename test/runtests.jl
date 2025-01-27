import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Muninn
using Test
using JLD2
using Infiltrator
using Dates

include("TI.jl")
include("MB.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Construct TI models by default" TI_creation_default_test()
@testset "Construct TI models with input values" TI_creation_values_test()
@testset "Apply MB model" apply_MB_test()
