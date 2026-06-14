using Sleipnir, Muninn
using Statistics, Printf

RGI = "RGI60-11.01450"
TSPAN = (2000.0, 2020.0)

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI => rgi_paths[RGI]),
    use_MB = true, use_velocities = false, climate_data_source = :W5E5))
glacier = Sleipnir.initialize_glaciers([RGI], params)[1]

# Check August 2010
step = Sleipnir.Float(1.0/12.0)
t_aug = Sleipnir.Float(2010.0 + 8/12)
Sleipnir.get_cumulative_climate!(glacier.climate, t_aug, step)
cs = glacier.climate.climate_step
crs = glacier.climate.climate_raw_step

temp_raw = vec(crs.temp.data)
prcp_raw = vec(crs.prcp.data)
grad_raw = vec(crs.gradient.data)

@printf("August 2010 raw climate:\n")
@printf("  ref_hgt: %.0f m\n", cs.ref_hgt)
@printf("  Days: %d\n", length(temp_raw))
@printf("  Temp: %.2f to %.2f C (mean: %.2f)\n",
    minimum(temp_raw), maximum(temp_raw), mean(temp_raw))
@printf("  Grad: %.6f to %.6f C/m (mean: %.6f)\n",
    minimum(grad_raw), maximum(grad_raw), mean(grad_raw))

# Manual computation at 2500m
h_test = 2500.0
ΔS = h_test - cs.ref_hgt

pdd_calc = 0.0
snow_calc = 0.0
for d in eachindex(temp_raw)
    T_d = temp_raw[d]
    g_d = grad_raw[d]
    p_d = prcp_raw[d]
    T_eff = T_d  # temp_bias = 0
    pdd_day = max(0, T_eff + g_d * ΔS)
    snow_day = p_d * clamp((2 - T_eff - g_d * ΔS) / 2, 0, 1)
    pdd_calc += pdd_day
    snow_calc += snow_day
end

@printf("\nManual computation at %.0f m:\n", h_test)
@printf("  ΔS: %.0f m\n", ΔS)
@printf("  PDD: %.2f C-days\n", pdd_calc)
@printf("  Snow: %.2f mm\n", snow_calc)

# Now check what downscale_2D_climate gives
cd_2D = Sleipnir.downscale_2D_climate(
    glacier.climate.climate_step,
    glacier.climate.climate_raw_step,
    glacier.S, glacier.Coords; include_topography = false)

# Find the grid point closest to h_test
elev_diff = glacier.S .- cs.ref_hgt
idx = argmin(abs.(elev_diff .- (h_test - cs.ref_hgt)))

@printf("\nDownscaled 2D at nearest grid point:\n")
@printf("  Actual elevation: %.0f m\n", glacier.S[idx])
@printf("  Actual ΔS: %.0f m\n", glacier.S[idx] - cs.ref_hgt)
@printf("  PDD: %.2f C-days\n", cd_2D.PDD[idx])
@printf("  Snow: %.2f mm\n", cd_2D.snow[idx])

# Compare
@printf("\nComparison:\n")
@printf("  PDD ratio: %.2f%%\n", cd_2D.PDD[idx] / pdd_calc * 100)
@printf("  Snow ratio: %.2f%%\n", cd_2D.snow[idx] / snow_calc * 100)
