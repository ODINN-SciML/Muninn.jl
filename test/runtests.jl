import Pkg
function is_included_in_repl()
    for frame in StackTraces.stacktrace()
        if occursin("start_repl_backend", string(frame.func))
            return true
        end
    end
    return false
end

Pkg.activate(dirname(Base.current_project()))
Pkg.instantiate() # Need this to setup the ODINN env for multiprocessing
if is_included_in_repl()
    # The Project.toml of the test environment to be used when running with include is in a subfolder to avoid that Julia uses this file in test mode
    Pkg.activate(dirname(Base.current_project())*"/test/test_env/")
    Pkg.resolve()
end

if !parse(Bool, get(ENV, "CI", "false"))
    using Revise
end
using Muninn
using Sleipnir: Parameters, Model
using Test
using JLD2
using Infiltrator
using Dates
using JET
using Aqua

include("TI.jl")
include("MB.jl")
include("Aqua.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Run all tests" begin
    @testset "Construct TI models by default" TI_creation_default_test()
    @testset "Construct TI models with input values" TI_creation_values_test()
    @testset "Synthetic TI MB field" TI_synthetic_field_test()
    @testset "MB compatibility helpers" model_compatibility_utils_test()
    @testset "Apply MB model" apply_MB_test()
    @testset "Aqua" test_Aqua()
end
