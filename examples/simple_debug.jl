using Sleipnir, Muninn
using Statistics, Dates, Printf

RGI = "RGI60-11.01450"
TSPAN = (2000.0, 2020.0)

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI => rgi_paths[RGI]),
    use_MB = true, use_velocities = false, climate_data_source = :W5E5))
glacier = Sleipnir.initialize_glaciers([RGI], params)[1]

ice = .!glacier.mask
@printf("Glacier: %s, cells: %d, elev: %.0f-%.0f m\n",
    glacier.rgi_id, count(ice), minimum(glacier.S[ice]), maximum(glacier.S[ice]))

# Check August 2010
step = Sleipnir.Float(1.0/12.0)
t_aug = Sleipnir.Float(2010.0 + 8/12)
Sleipnir.get_cumulative_climate!(glacier.climate, t_aug, step)
cs = glacier.climate.climate_step
crs = glacier.climate.climate_raw_step

@printf("\nAugust 2010:\n")
@printf("  avg_temp at ref_hgt: %.2f C\n", cs.avg_temp)
@printf("  total PDD: %.2f C-days\n", cs.temp)
@printf("  total prcp: %.2f mm\n", cs.prcp)

# Downscale
cd_2D = Sleipnir.downscale_2D_climate(
    glacier.climate.climate_step,
    glacier.climate.climate_raw_step,
    glacier.S, glacier.Coords; include_topography = false)

@printf("\nDownscaled 2D for August:\n")
@printf("  PDD mean on glacier: %.2f C-days\n", mean(cd_2D.PDD[ice]))
@printf("  Snow mean on glacier: %.2f mm\n", mean(cd_2D.snow[ice]))

# Compute MB with defaults
model_def = TImodel1(params; DDF = Sleipnir.Float(0.007), prcp_fac = Sleipnir.Float(1.0))
t0, t1 = Sleipnir.Float(TSPAN[1]), Sleipnir.Float(TSPAN[2])
mb_default = Muninn.compute_mean_annual_MB(model_def, glacier, t0, t1)
@printf("\nMB with DDF=0.007, prcp_fac=1.0: %.3f m/yr\n", mb_default)
