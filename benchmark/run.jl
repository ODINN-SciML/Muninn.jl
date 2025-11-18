import Pkg
Pkg.activate(dirname(Base.current_project()))

using Muninn
using BenchmarkTools
using Logging
Logging.disable_logging(Logging.Info)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10

println("# Performance benchmark")

rgi_ids = ["RGI60-11.03638"]
rgi_paths = get_rgi_paths()
params = Sleipnir.Parameters(
    simulation = SimulationParameters(
    use_MB = true,
    use_velocities = false,
    tspan = (2010.0, 2015.0),
    test_mode = true,
    rgi_paths = rgi_paths),
)
glaciers = initialize_glaciers(rgi_ids, params)
TI1 = TImodel1(params)
model = Sleipnir.Model(nothing, TI1, nothing) # This test only needs a mass balance model
glacier = initialize_glaciers(rgi_ids, params)[1]
t = 2015.0
step_MB = 1.0/12.0

println("## Benchmark of MB_timestep")
t = @benchmark MB_timestep($model, $glacier, $step_MB, $t)
display(t)
println("")
