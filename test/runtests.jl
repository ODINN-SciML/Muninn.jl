import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Test
using JLD2
using Infiltrator

using Huginn

# include("PDE_UDE_solve.jl")
include("halfar.jl")
# include("mass_conservation.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

atol = 0.01
# @testset "PDE and UDE SIA solvers without MB" pde_solve_test(atol; MB=false, fast=true)

atol = 2.0
# @testset "PDE and UDE SIA solvers with MB" pde_solve_test(atol; MB=true, fast=true)

@testset "Halfar Solution" halfar_test(rtol=0.02, atol=1.0)

# @testset "Conservation of Mass - Flat Bed" unit_mass_flatbed_test(rtol=1.0e-7, atol=1000)