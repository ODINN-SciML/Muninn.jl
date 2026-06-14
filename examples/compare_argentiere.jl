#=
Argentière (RGI60-11.03638) — Muninn calibration + per-elevation MB profile
=============================================================================
Calibrates TImodel1 (prcp_fac=2.5, t_melt=0, W5E5) against Hugonnet geodetic
MB, then computes per-elevation-band annual MB for comparison with OGGM output.
=#

using Muninn, Sleipnir
using Statistics, Printf

const RGI = "RGI60-11.03638"
const TSPAN = (2000.0, 2020.0)
const NYEARS = Float64(TSPAN[2] - TSPAN[1])

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI => rgi_paths[RGI]),
    use_MB = true, use_velocities = false, climate_data_source = :W5E5))

glacier = Sleipnir.initialize_glaciers([RGI], params)[1]

println("="^70)
println("Muninn TImodel1 — Argentière ($RGI)")
println("="^70)
@printf("  ref_hgt=%.0fm   ice cells=%d\n",
    glacier.climate.ref_hgt, count(.!glacier.mask))
@printf("  t_melt=0°C (Muninn default)   prcp_fac=2.5 (Step 1 fixed)\n")
@printf("  Geodetic target: %.4f m w.e. yr⁻¹\n", glacier.geodetic_MB)
println()

# ── Calibrate ─────────────────────────────────────────────────────────────────
cal = calibrate_ti_model(glacier, params)
t0, t1 = Sleipnir.Float(TSPAN[1]), Sleipnir.Float(TSPAN[2])
mb_cal = compute_mean_annual_MB(cal, glacier, t0, t1)

@printf("  Calibrated DDF   : %.4f mm w.e. °C⁻¹ d⁻¹\n", cal.DDF * 1000)
@printf("  prcp_fac         : %.4f\n", cal.prcp_fac)
@printf("  Glacier-wide MB  : %.4f m w.e. yr⁻¹  (error=%.2e)\n",
    mb_cal, mb_cal - glacier.geodetic_MB)
println()

# ── Per-elevation MB profile using calibrated model ───────────────────────────
total_mb = Muninn.compute_cumulative_MB(cal, glacier, t0, t1)
ann_mb = total_mb ./ NYEARS   # m w.e. / yr per cell

ice = .!glacier.mask
S_ice = glacier.S[ice]
mb_ice = ann_mb[ice]

BIN = 50
lo_edge = floor(Int, minimum(S_ice) / BIN) * BIN
hi_edge = ceil(Int, maximum(S_ice) / BIN) * BIN

# Accumulation and melt separately
step = Sleipnir.Float(1.0/12.0)
total_acc = zeros(Sleipnir.Float, size(glacier.S))
total_melt = zeros(Sleipnir.Float, size(glacier.S))
t = t0 + step
while t <= t1 + step/2
    Sleipnir.get_cumulative_climate!(glacier.climate, t, step)
    cd = Sleipnir.downscale_2D_climate(
        glacier.climate.climate_step, glacier.climate.climate_raw_step,
        glacier.S, glacier.Coords; include_topography = false,
        temp_bias = get_temp_bias(cal))
    total_acc .+= 0.001f0 .* cal.prcp_fac .* cd.snow
    total_melt .+= cal.DDF .* cd.PDD
    global t += step
end
acc_ice = (total_acc ./ NYEARS)[ice]
melt_ice = (total_melt ./ NYEARS)[ice]

@printf("  %-18s %8s %12s %10s %10s\n",
    "Elev band (m)", "n cells", "MB (m/yr)", "accum", "melt")
println("  " * "-"^62)
for lo in lo_edge:BIN:(hi_edge - BIN)
    hi = lo + BIN
    idx = findall(s -> lo <= s < hi, S_ice)
    isempty(idx) && continue
    @printf("  %4d–%-5d       %8d %12.3f %10.3f %10.3f\n",
        lo, hi, length(idx),
        mean(mb_ice[idx]), mean(acc_ice[idx]), mean(melt_ice[idx]))
end
