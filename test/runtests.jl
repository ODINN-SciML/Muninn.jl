import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Test
using JLD2
using Infiltrator

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"
