#=
Debug Gulkana Step-1 calibration: why does it fall back to Step 2?
=#

using Muninn, Sleipnir, Statistics, Printf

RGI = "RGI60-01.00570"
TSPAN = (2000.0, 2020.0)

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI=>rgi_paths[RGI]), use_MB = true, use_velocities = false,
    climate_data_source = :W5E5))

println("Initializing Gulkana ($RGI)...")
glacier = Sleipnir.initialize_glaciers([RGI], params)[1]

# Compute winter prcp_fac
prcp_fac_winter = Sleipnir.get_winter_prcp_factor(glacier, params)
println("Winter-derived prcp_fac: $(round(prcp_fac_winter; digits=4))")

# Target MB
mb_target = Sleipnir.Float(glacier.geodetic_MB)
println("Geodetic MB target: $(round(mb_target; digits=4)) m w.e. yr⁻¹\n")

# Step-1 setup
DDF_min = Sleipnir.Float(params.physical.DDF_min)
DDF_max = Sleipnir.Float(params.physical.DDF_max)
tb_zero = Sleipnir.Float(0.0)
t_start = Sleipnir.Float(TSPAN[1])
t_end = Sleipnir.Float(TSPAN[2])
step = Sleipnir.Float(1.0/12.0)

println("Step 1: Bracketing check with winter prcp_fac = $(round(prcp_fac_winter; digits=4))")
println("  DDF search interval: [$(round(DDF_min*1000; digits=2)), $(round(DDF_max*1000; digits=2))] mm/°C/d")

# Residuals at bounds
function residual_ddf(DDF)
    model_trial = Muninn.TImodel1{Sleipnir.Float}(DDF, prcp_fac_winter, tb_zero)
    return Muninn.compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) -
           mb_target
end

r_min = residual_ddf(DDF_min)
r_max = residual_ddf(DDF_max)

@printf("  r(DDF_min=%.2f) = %.4f m w.e./yr\n", DDF_min*1000, r_min)
@printf("  r(DDF_max=%.2f) = %.4f m w.e./yr\n", DDF_max*1000, r_max)

if r_min * r_max <= 0
    println("  ✓ Bracketing SUCCESSFUL — Step 1 solves\n")
else
    println("  ✗ Bracketing FAILED — Step 1 would fall back to Step 2\n")

    # Diagnose which direction the residuals are
    if r_min > 0 && r_max > 0
        println("  DIAGNOSIS: Both residuals positive (both too negative MB)")
        println("    → All DDF values give net accumulation; target is too negative")
        println("    → Need lower prcp_fac or higher temp_bias")
    elseif r_min < 0 && r_max < 0
        println("  DIAGNOSIS: Both residuals negative (both too positive MB)")
        println("    → All DDF values give net ablation; target is too positive")
        println("    → Need higher prcp_fac or lower temp_bias")
    end
end

# Compare with fixed 2.5
println("\n" * "="^60)
println("Comparison: if Step 1 used fixed prcp_fac = 2.5 instead")
println("="^60)

prcp_fac_fixed = Sleipnir.Float(2.5)
function residual_ddf_fixed(DDF)
    model_trial = Muninn.TImodel1{Sleipnir.Float}(DDF, prcp_fac_fixed, tb_zero)
    return Muninn.compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) -
           mb_target
end

r_min_fixed = residual_ddf_fixed(DDF_min)
r_max_fixed = residual_ddf_fixed(DDF_max)

@printf("  r(DDF_min=%.2f) = %.4f m w.e./yr  (was %.4f)\n", DDF_min*1000, r_min_fixed,
    r_min)
@printf("  r(DDF_max=%.2f) = %.4f m w.e./yr  (was %.4f)\n", DDF_max*1000, r_max_fixed,
    r_max)

if r_min_fixed * r_max_fixed <= 0
    println("  ✓ Fixed 2.5 would BRACKET\n")
else
    println("  ✗ Fixed 2.5 also wouldn't bracket\n")
end

# Show specific MB at glacier-mean elevation (climate-based, ignoring topography)
println("\nMB at mean glacier elevation (climate_step summary):")
Sleipnir.get_cumulative_climate!(glacier.climate, Sleipnir.Float(2010.0), Sleipnir.Float(1.0))
cs = glacier.climate.climate_step
@printf("  Mean annual precip:  %.1f mm w.e.\n", cs.prcp)
@printf("  Mean annual PDDs:    %.1f °C·days\n", cs.temp)
@printf("  Avg temp gradient:   %.4f °C/m\n", cs.avg_gradient)
