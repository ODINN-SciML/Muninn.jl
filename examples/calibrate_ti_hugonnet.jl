#=
Calibrate Temperature-Index (TI) Models with Hugonnet Geodetic Mass Balance
=========================================================================

This script demonstrates how to calibrate TI models (TImodel1 with 1 DDF)
for WGMS-monitored glaciers using geodetic mass balance observations from Hugonnet et al. 2021.

Structure follows the ODINN forward simulation tutorial pattern.

Glaciers selected (7 WGMS reference glaciers, globally distributed):
  - RGI60-11.03638: Argentière (13.8 km²) - Mont Blanc, Alps
  - RGI60-11.01450: Aletschgletscher (82.2 km²) - Valais, Alps
  - RGI60-11.01238: Rhonegletscher (15.8 km²) - Central Alps
  - RGI60-11.00897: Hintereisferner (8.0 km²) - Ötztal, Austria
  - RGI60-08.00213: Storglaciaeren (3.4 km²) - Kebnekaise, Sweden
  - RGI60-01.00570: Gulkana (17.6 km²) - Alaska Range, USA
  - RGI60-02.05098: Peyto (9.7 km²) - Canadian Rockies

Hugonnet et al. 2021 provides glacier-wide geodetic mass balance for 2000-2020.
The calibration adjusts the degree-day factor (DDF) to match observed MB,
using a glacier-specific precipitation factor derived from winter precipitation.

Usage:
  julia calibrate_ti_hugonnet.jl

Output:
  - Prints calibrated DDF and prcp_fac for each glacier
  - Calibrated models can be used for forward simulations
=#

using Revise
using Muninn
using Sleipnir
using Printf

# ============================================================================
# Step 1: Parameter Initialization (follows forward simulation tutorial)
# ============================================================================

# WGMS reference glaciers — globally distributed, with continuous monitoring
const WGMS_GLACIERS = [
    "RGI60-11.03638",  # Argentière (France)
    "RGI60-11.01450",  # Aletschgletscher (Switzerland)
    "RGI60-11.01238",  # Rhonegletscher (Switzerland)
    "RGI60-11.00897",  # Hintereisferner (Austria)
    "RGI60-08.00213",  # Storglaciaeren (Sweden)
    "RGI60-01.00570",  # Gulkana (Alaska)
    "RGI60-02.05098"  # Peyto (Canada)
]

# Hugonnet observation period: 2000-01-01 to 2020-01-01
const CALIBRATION_TSPAN = (2000.0, 2020.0)

# Standalone Muninn runs single-process (multiprocessing is handled by Huginn/ODINN).
# At scale (hundreds/thousands of glaciers), use calibrate_MB_model via a
# Huginn Prediction or ODINN Inversion — their Parameters constructor loads Muninn
# on the workers so the pmap in calibrate_MB_model distributes automatically.
params = Sleipnir.Parameters(
    simulation = Sleipnir.SimulationParameters(
    tspan = CALIBRATION_TSPAN,
    multiprocessing = false,
    workers = 1,
    rgi_paths = Sleipnir.get_rgi_paths(),
    use_MB = true,
    use_velocities = false,
    test_mode = false
)
)

# ============================================================================
# Step 2: Glacier Initialization
# ============================================================================

@info "Initializing glaciers with Hugonnet geodetic MB data..."

# Initialize glaciers - this automatically loads Hugonnet dhdtData
# via Sleipnir's _default_hugonnet_dhdt() function
# The glacier.dhdtData field contains the observation period and dh/dt
# The glacier.geodetic_MB field contains the glacier-wide mean MB in m w.e. yr⁻¹
glaciers = initialize_glaciers(WGMS_GLACIERS, params)

# Verify all glaciers have geodetic MB data
for glacier in glaciers
    if isnothing(glacier.dhdtData) || !isfinite(glacier.geodetic_MB)
        @warn "Glacier $(glacier.rgi_id) is missing Hugonnet geodetic MB data. Skipping."
    end
end

@info "All glaciers initialized successfully with Hugonnet data."

# ============================================================================
# Step 3: TI Model Calibration
# ============================================================================

@info "Calibrating TI models with 1 DDF for each glacier..."

# calibrate_ti_model() performs a 3-step cascade for a single glacier:
# 1. Calibrate DDF with prcp_fac derived from winter precipitation (default: :from_winter_prcp)
# 2. If DDF can't bracket target, calibrate prcp_fac with DDF at boundary
# 3. If both fail, calibrate temp_bias with DDF and prcp_fac at boundaries
#
# By default, prcp_fac is glacier-specific and derived from mean winter daily precipitation,
# following OGGM's approach. Pass prcp_fac=2.5 to use a fixed global value instead.
#
# It returns a fresh TImodel1 with calibrated (DDF, prcp_fac, temp_bias).
# To calibrate the MB model of a whole simulation, use
# calibrate_MB_model(model, glaciers, params) instead (or just set
# `calibrate_MB = true` in SimulationParameters).

# TI1 calibration uses static glacier geometry, so no iceflow model is needed
# (the Model.iceflow field is never read during calibration).
# calibrate_MB_model expands the single TImodel1 template into a per-glacier
# vector of calibrated models, running the fits under pmap — distributed across
# workers when multiprocessing = true.
model = Model(; iceflow = nothing, mass_balance = TImodel1(params))
model = calibrate_MB_model(model, glaciers, params)
calibrated_models = model.mass_balance

# ============================================================================
# Step 4: Results Display and Validation
# ============================================================================

println(repeat("=", 80))
println("TI Model Calibration Results (Hugonnet 2000-2020)")
println(repeat("=", 80))

# Create a default model for comparison
default_model = TImodel1(params)

for (i, (glacier, cal_model)) in enumerate(zip(glaciers, calibrated_models))
    println("\n[$(i)] Glacier: $(glacier.rgi_id)")

    # Display input data
    println("  📊 Input:")
    glacier_cells = count(.!glacier.mask)
    area_m2 = glacier_cells * abs(glacier.Δx * glacier.Δy)
    println("    Area: $(round(area_m2 / 1e6; digits=2)) km²")
    println("    Hugonnet geodetic MB: $(round(glacier.geodetic_MB; digits=4)) m w.e. yr⁻¹")

    # Display default model MB for comparison
    default_mb = compute_mean_annual_MB(
        default_model, glacier,
        Sleipnir.Float(CALIBRATION_TSPAN[1]),
        Sleipnir.Float(CALIBRATION_TSPAN[2])
    )
    println("    Default model MB: $(round(default_mb; digits=4)) m w.e. yr⁻¹")

    # Display calibrated parameters
    println("  ✅ Calibrated TImodel1:")
    println("    DDF: $(round(cal_model.DDF * 1000; digits=4)) mm w.e. °C⁻¹ d⁻¹")
    println("    prcp_fac: $(cal_model.prcp_fac)")
    println("    temp_bias: $(round(cal_model.temp_bias; digits=4)) °C")

    # Validate: compute MB with calibrated model
    cal_mb = compute_mean_annual_MB(
        cal_model, glacier,
        Sleipnir.Float(CALIBRATION_TSPAN[1]),
        Sleipnir.Float(CALIBRATION_TSPAN[2])
    )
    mb_error = cal_mb - glacier.geodetic_MB
    println("  🎯 Validation:")
    println("    Calibrated model MB: $(round(cal_mb; digits=4)) m w.e. yr⁻¹")
    println("    Error vs Hugonnet: $(round(mb_error; digits=6)) m w.e. yr⁻¹")
    println("    Relative error: $(round(abs(mb_error / glacier.geodetic_MB) * 100; digits=2))%")
end

println(repeat("=", 80))
println("Calibration complete.")
println("Models are ready for use in forward simulations.")
println(repeat("=", 80))

# ============================================================================
# Step 5: Cumulative MB maps
# ============================================================================

@info "Plotting cumulative mass balance maps..."

t_start = Sleipnir.Float(CALIBRATION_TSPAN[1])
t_end = Sleipnir.Float(CALIBRATION_TSPAN[2])

for (glacier, cal_model) in zip(glaciers, calibrated_models)
    fig = plot_cumulative_mb(cal_model, glacier, t_start, t_end;
        title = "Mean Annual MB $(Int(CALIBRATION_TSPAN[1]))–$(Int(CALIBRATION_TSPAN[2]))")
    path = save_figure(fig, joinpath("plots", "cumulative_mb_$(glacier.rgi_id).png"))
    @info "Saved: $path"
end

# ============================================================================
# Step 6: Save Results (Optional)
# ============================================================================

# Uncomment to save calibrated models to a JLD2 file
# using JLD2
# @save "ti_calibrated_hugonnet.jld2" calibrated_models WGMS_GLACIERS CALIBRATION_TSPAN
