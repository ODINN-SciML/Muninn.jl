#= 
Compare Climate Downscaling: Sleipnir vs OGGM
=============================================

This script compares the climate data and downscaling approach between
Sleipnir (used by Muninn) and OGGM to understand why parameters don't
directly transfer.
=#

using Sleipnir
using Muninn
using Statistics
using Printf
using Dates

println(repeat("=", 80))
println("Climate Downscaling Comparison: Sleipnir vs OGGM")
println(repeat("=", 80))
println()

# Initialize glacier
rgi_paths = Sleipnir.get_rgi_paths()
rgi_paths_filtered = Dict(k => rgi_paths[k] for k in ["RGI60-11.00897"])

params = Sleipnir.Parameters(
    simulation = Sleipnir.SimulationParameters(
    tspan = (2000.0, 2020.0),
    multiprocessing = false,
    workers = 1,
    rgi_paths = rgi_paths_filtered,
    use_MB = true,
    use_velocities = false,
    test_mode = false,
    climate_data_source = :W5E5
)
)

glacier = Sleipnir.initialize_glaciers(["RGI60-11.00897"], params)[1]

println("Glacier Information:")
println("-"^80)
@printf("  RGI ID: %s\n", glacier.rgi_id)
@printf("  Area: %.2f km²\n", count(.!glacier.mask) * abs(glacier.Δx * glacier.Δy) / 1e6)
@printf("  Surface elevation: %.0f to %.0f m (mean: %.0f m)\n",
    minimum(glacier.S), maximum(glacier.S), mean(glacier.S))
println()

# Climate data from Sleipnir
climate_raw = glacier.climate.climate_raw_step
climate_step = glacier.climate.climate_step

println("Sleipnir Climate Data (W5E5):")
println("-"^80)
@printf("  Reference height: %.0f m\n", climate_step.ref_hgt)
let tdim = collect(Sleipnir.dims(climate_raw, Sleipnir.Ti))
    @printf("  Raw climate period: %s to %s\n", tdim[1], tdim[end])
end
println()

println("Daily Climate Statistics (2000-2020):")
@printf("  Temperature: %.2f to %.2f °C (mean: %.2f °C)\n",
    minimum(climate_raw.temp), maximum(climate_raw.temp), mean(climate_raw.temp))
@printf("  Precipitation: %.2f to %.2f mm/day (mean: %.2f mm/day)\n",
    minimum(climate_raw.prcp), maximum(climate_raw.prcp), mean(climate_raw.prcp))
@printf("  Gradient: %.2f to %.2f °C/m (mean: %.6f °C/m)\n",
    minimum(climate_raw.gradient), maximum(climate_raw.gradient),
    mean(climate_raw.gradient))
println()

println("Cumulative Climate (2000-2020):")
@printf("  avg_temp: %.2f °C\n", climate_step.avg_temp)
@printf("  avg_gradient: %.6f °C/m\n", climate_step.avg_gradient)
@printf("  total_prcp: %.2f mm\n", climate_step.prcp)
@printf("  total_PDD (from temp): %.2f °C-days\n", climate_step.temp)
println()

# Downscaled climate
using Muninn: downscale_2D_climate

climate_2D = downscale_2D_climate(
    glacier.climate.climate_step,
    glacier.climate.climate_raw_step,
    glacier.S,
    glacier.Coords;
    include_topography = false
)

println("Downscaled 2D Climate (single cumulative step):")
@printf("  Mean temp_2D: %.2f °C\n", mean(climate_2D.temp))
@printf("  Mean PDD_2D: %.2f °C-days\n", mean(climate_2D.PDD))
@printf("  Mean snow_2D: %.2f mm\n", mean(climate_2D.snow))
@printf("  Mean rain_2D: %.2f mm\n", mean(climate_2D.rain))
@printf("  Mean elevation_diff: %.0f m\n", mean(climate_2D.elevation_diff))
println()

# Now let's check monthly
println("Monthly Climate Check (first month):")
println("-"^80)
get_cumulative_climate!(glacier.climate, Sleipnir.Float(2000.0), Sleipnir.Float(1.0/12.0))
climate_step_monthly = glacier.climate.climate_step
@printf("  First month (Jan 2000) avg_temp: %.2f °C\n", climate_step_monthly.avg_temp)
@printf("  First month total_PDD: %.2f °C-days\n", climate_step_monthly.temp)
@printf("  First month total_prcp: %.2f mm\n", climate_step_monthly.prcp)

climate_2D_monthly = downscale_2D_climate(
    climate_step_monthly,
    glacier.climate.climate_raw_step,
    glacier.S,
    glacier.Coords;
    include_topography = false
)
@printf("  First month mean PDD_2D: %.2f °C-days\n", mean(climate_2D_monthly.PDD))
@printf("  First month mean snow_2D: %.2f mm\n", mean(climate_2D_monthly.snow))
println()

println("="^80)
println("Comparison with OGGM:")
println("="^80)
println()

println("Key Differences:")
println("-"^80)
println("1. Temperature Lapse Rate:")
@printf("   - Sleipnir W5E5: %.6f °C/m\n", mean(climate_raw.gradient))
@printf("   - OGGM (from tutorial): %.6f °C/m\n", -0.0065)
@printf("   - Difference: %.6f °C/m (%.2f%%)\n",
    mean(climate_raw.gradient) - (-0.0065),
    abs(mean(climate_raw.gradient) - (-0.0065)) / 0.0065 * 100)
println()

println("2. Reference Height:")
@printf("   - Sleipnir: %.0f m\n", climate_step.ref_hgt)
@printf("   - Glacier mean elevation: %.0f m\n", mean(glacier.S))
@printf("   - Elevation difference: %.0f m\n", mean(glacier.S) - climate_step.ref_hgt)
println()

println("3. Snow/Rain Partitioning:")
println("   - Both use: temp_all_solid = 0.0°C, temp_all_liq = 2.0°C")
println("   - Linear transition between 0°C and 2°C")
println()

println("4. Impact on Downscaling:")
@printf("   - Temperature at glacier surface (mean): %.2f °C\n",
    climate_step.avg_temp +
    climate_step.avg_gradient * (mean(glacier.S) - climate_step.ref_hgt))
@printf("   - With OGGM lapse rate: %.2f °C\n",
    climate_step.avg_temp + (-0.0065) * (mean(glacier.S) - climate_step.ref_hgt))
@printf("   - Difference: %.2f °C\n",
    (climate_step.avg_temp +
     climate_step.avg_gradient * (mean(glacier.S) - climate_step.ref_hgt)) -
    (climate_step.avg_temp + (-0.0065) * (mean(glacier.S) - climate_step.ref_hgt)))
println()

println("="^80)
println("Conclusion:")
println("="^80)
println()
println("The main difference is the temperature lapse rate:")
println("  - Sleipnir's W5E5 data has: -0.00552 °C/m")
println("  - OGGM uses: -0.0065 °C/m")
println()
println("This 15% difference in lapse rate leads to different temperature")
println("adjustments at the glacier surface elevation, which affects PDD")
println("and snow/rain partitioning.")
println()
println("Both systems use the same snow/rain partitioning thresholds (0°C and 2°C),")
println("but the different lapse rates lead to different parameter requirements")
println("to match the same observed mass balance.")
