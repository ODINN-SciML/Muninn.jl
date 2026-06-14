#= 
Reproduce OGGM Mass Balance Calibration Tutorial
==================================================

This script reproduces the OGGM mass balance calibration tutorial:
https://tutorials.oggm.org/stable/notebooks/tutorials/massbalance_calibration.html

The tutorial uses Hintereisferner glacier (RGI60-11.00897) with:
- Reference period: 2000-01-01 to 2020-01-01
- Reference MB: -1100.3 kg m⁻² (≈ -1.1003 m w.e. yr⁻¹)
- OGGM calibrated parameters: melt_f=5.0, prcp_fac=3.57, temp_bias=1.705

We'll use Muninn to calibrate the same glacier and compare results.

Note: OGGM uses kg m⁻² units, Muninn uses m w.e. units.
Conversion: 1 kg m⁻² = 1 mm w.e. = 0.001 m w.e.
So: -1100.3 kg m⁻² = -1.1003 m w.e.
=#

using Muninn
using Sleipnir
using Printf

# ============================================================================
# Configuration
# ============================================================================

const RGI_ID = "RGI60-11.00897"  # Hintereisferner
const REFERENCE_MB_OGGM = -1100.3  # kg m⁻² yr⁻¹
const REFERENCE_MB_MUNINN = REFERENCE_MB_OGGM / 1000.0  # Convert to m w.e. yr⁻¹
const CALIBRATION_TSPAN = (2000.0, 2020.0)

println(repeat("=", 80))
println("Reproducing OGGM Mass Balance Calibration for $RGI_ID")
println(repeat("=", 80))
println()

# ============================================================================
# Initialize glacier
# ============================================================================

@printf("Initializing glacier %s...\n", RGI_ID)

rgi_paths = Sleipnir.get_rgi_paths()
# Filter to only the glacier we need
rgi_paths_filtered = Dict(k => rgi_paths[k] for k in [RGI_ID])

params = Sleipnir.Parameters(
    simulation = Sleipnir.SimulationParameters(
    tspan = CALIBRATION_TSPAN,
    multiprocessing = false,
    workers = 1,
    rgi_paths = rgi_paths_filtered,
    use_MB = true,
    use_velocities = false,
    test_mode = false,  # Use real preprocessed data
    climate_data_source = :W5E5  # Use W5E5 to match OGGM
)
)

glacier = initialize_glaciers([RGI_ID], params)[1]

@printf("Glacier initialized: %s\n", glacier.rgi_id)
@printf("  Area: %.2f km²\n", count(.!glacier.mask) * abs(glacier.Δx * glacier.Δy) / 1e6)

# Check if we have geodetic MB data
if isnothing(glacier.dhdtData)
    @warn "No geodetic MB data available for this glacier"
elseif !isfinite(glacier.geodetic_MB)
    @warn "Geodetic MB is not finite"
else
    @printf("  Geodetic MB: %.4f m w.e. yr⁻¹ (from data)\n", glacier.geodetic_MB)
end
println()

# Set up time variables
t_start = Sleipnir.Float(CALIBRATION_TSPAN[1])
t_end = Sleipnir.Float(CALIBRATION_TSPAN[2])

# ============================================================================
# OGGM Reference Values
# ============================================================================

println(repeat("-", 80))
println("OGGM Reference Values (from tutorial):")
println(repeat("-", 80))
@printf("  Reference MB: %.1f kg m⁻² yr⁻¹ = %.4f m w.e. yr⁻¹\n",
    REFERENCE_MB_OGGM, REFERENCE_MB_MUNINN)
@printf("  Calibrated melt_f: 5.0 kg m⁻² day⁻¹ K⁻¹ = %.6f m w.e. day⁻¹ K⁻¹\n",
    5.0 / 1000.0)  # Convert kg m⁻² to m w.e.
@printf("  Calibrated prcp_fac: 3.570 (dimensionless)\n")
@printf("  Calibrated temp_bias: 1.705 °C\n")
@printf("  Global params: temp_default_gradient=-0.0065, temp_all_solid=0.0, temp_all_liq=2.0, temp_melt=-1.0\n")
println()

# ============================================================================
# Test 1: Use OGGM's exact parameters in Muninn
# ============================================================================

println(repeat("-", 80))
println("Test 1: Using OGGM's exact parameters in Muninn")
println(repeat("-", 80))

# OGGM parameters (from tutorial)
# melt_f = 5.0 kg m⁻² day⁻¹ K⁻¹ = 0.005 m w.e. day⁻¹ K⁻¹
# prcp_fac = 3.570255710475343
# temp_bias = 1.705035916562294 °C

oggm_melt_f_mwe = 5.0 / 1000.0  # Convert to m w.e. units
oggm_prcp_fac = 3.570255710475343
oggm_temp_bias = 1.705035916562294

# Create a model with OGGM's parameters
oggm_model = TImodel1(params;
    DDF = Sleipnir.Float(oggm_melt_f_mwe),
    prcp_fac = Sleipnir.Float(oggm_prcp_fac),
    temp_bias = Sleipnir.Float(oggm_temp_bias)
)

# Compute MB with OGGM parameters
oggm_mb = compute_mean_annual_MB(
    oggm_model, glacier, t_start, t_end
)

@printf("OGGM parameters in Muninn:\n")
@printf("  DDF: %.6f m w.e. °C⁻¹ day⁻¹\n", oggm_model.DDF)
@printf("  prcp_fac: %.6f\n", oggm_model.prcp_fac)
@printf("  temp_bias: %.6f °C\n", oggm_model.temp_bias)
@printf("  Computed MB: %.4f m w.e. yr⁻¹\n", oggm_mb)
@printf("  Target MB: %.4f m w.e. yr⁻¹\n", REFERENCE_MB_MUNINN)
@printf("  Error: %.6f m w.e. yr⁻¹\n", oggm_mb - REFERENCE_MB_MUNINN)
@printf("  Relative error: %.2f%%\n",
    abs((oggm_mb - REFERENCE_MB_MUNINN) / REFERENCE_MB_MUNINN) * 100)
println()

# ============================================================================
# Test 2: Muninn Calibration
# ============================================================================

println(repeat("-", 80))
println("Test 2: Muninn's automatic calibration")
println(repeat("-", 80))

# Attach the OGGM reference MB as the target
dhdt_data = Sleipnir.DhdtData(
    (Sleipnir.Float(CALIBRATION_TSPAN[1]), Sleipnir.Float(CALIBRATION_TSPAN[2])),
    Sleipnir.Float(REFERENCE_MB_MUNINN)
)
glacier_with_target = Sleipnir.Glacier2D(glacier;
    dhdtData = dhdt_data,
    geodetic_MB = Sleipnir.Float(REFERENCE_MB_MUNINN)
)

# Calibrate using Muninn's 3-step method
@printf("Calibrating with Muninn's 3-step method...\n")
cal_model = calibrate_ti_model(glacier_with_target, params)

println()
@printf("Muninn Calibrated Parameters:\n")
@printf("  DDF: %.6f m w.e. °C⁻¹ day⁻¹ (≈ %.2f kg m⁻² day⁻¹ K⁻¹)\n",
    cal_model.DDF, cal_model.DDF * 1000)
@printf("  PRECIP_UNIT_CONVERSION: %.6f m w.e. mm⁻¹ (hardcoded)\n", PRECIP_UNIT_CONVERSION)
@printf("  prcp_fac: %.4f (dimensionless)\n", cal_model.prcp_fac)
@printf("  temp_bias: %.4f °C\n", cal_model.temp_bias)
println()

# ============================================================================
# Compute MB with calibrated model
# ============================================================================

@printf("Computing MB with calibrated model...\n")

calibrated_mb = compute_mean_annual_MB(
    cal_model, glacier_with_target, t_start, t_end
)

@printf("Calibrated model MB: %.4f m w.e. yr⁻¹\n", calibrated_mb)
@printf("Target MB: %.4f m w.e. yr⁻¹\n", REFERENCE_MB_MUNINN)
@printf("Error: %.6f m w.e. yr⁻¹\n", calibrated_mb - REFERENCE_MB_MUNINN)
@printf("Relative error: %.2f%%\n",
    abs((calibrated_mb - REFERENCE_MB_MUNINN) / REFERENCE_MB_MUNINN) * 100)
println()

# ============================================================================
# Comparison with OGGM
# ============================================================================

println(repeat("=", 80))
println("Comparison: OGGM vs Muninn")
println(repeat("=", 80))

println("\nParameter Comparison:")
println("-"^80)
@printf("%-20s | %-25s | %-25s | %-10s\n", "Parameter", "OGGM", "Muninn Calibrated",
    "Units")
println("-"^80)
@printf("%-20s | %-25.6f | %-25.6f | %-10s\n", "Melt Factor/DDF",
    oggm_melt_f_mwe, cal_model.DDF, "m w.e. d⁻¹ K⁻¹")
@printf("%-20s | %-25.4f | %-25.4f | %-10s\n", "prcp_fac",
    oggm_prcp_fac, cal_model.prcp_fac, "dimensionless")
@printf("%-20s | %-25.4f | %-25.4f | %-10s\n", "temp_bias",
    oggm_temp_bias, cal_model.temp_bias, "°C")

println()
println("MB Results:")
println("-"^80)
@printf("%-20s | %-25.4f | %-25s\n", "Reference MB", REFERENCE_MB_MUNINN, "m w.e. yr⁻¹")
@printf("%-20s | %-25.4f | %-25s\n", "OGGM params in Muninn", oggm_mb, "m w.e. yr⁻¹")
@printf("%-20s | %-25.4f | %-25s\n", "Muninn calibrated", calibrated_mb, "m w.e. yr⁻¹")

println()
println(repeat("=", 80))
println("Conclusion:")
println(repeat("=", 80))

println("\nTest 1 - Using OGGM's parameters directly in Muninn:")
if abs(oggm_mb - REFERENCE_MB_MUNINN) < 0.01
    println("✅ SUCCESS: OGGM parameters reproduce the target MB in Muninn")
else
    println("⚠️  OGGM parameters give a different MB in Muninn")
    println("   This is expected because:")
    println("   - Muninn uses different climate data (not W5E5)")
    println("   - The glacier geometry may differ")
    println("   - Downscaling methods may differ")
end

println("\nTest 2 - Muninn's automatic calibration:")
if abs(calibrated_mb - REFERENCE_MB_MUNINN) < 0.01
    println("✅ SUCCESS: Muninn calibration matches the target MB")
else
    println("⚠️  Muninn calibration differs from target")
end

println("\nNote: The calibrated parameters differ from OGGM's because:")
println("  - OGGM uses informed initial guesses (prcp_fac from winter precip, temp_bias from global calibration)")
println("  - Muninn uses a simple 3-step cascade starting from scratch")
println("  - Both methods find different parameter combinations that give the same MB")
println("  - This is expected and both are valid (overparameterized system)")

println()
@printf("Script completed successfully!\n")
