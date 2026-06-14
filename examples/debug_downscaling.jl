#=
Temperature reading debug — pinpoint where values diverge
=========================================================
Python ground truth for Aletsch raw_climate_(2000.0, 2020.0).nc:
  File covers:  1999-01-01 to 2019-12-31  (7670 days)
  temp mean=+3.9977°C  min=-17.68  max=+20.76
  Monthly means 2010: Jan=-6.3 Feb=-4.6 Mar=-2.1 Apr=+2.6 May=+4.6
                      Jun=+9.8 Jul=+13.3 Aug=+10.3 Sep=+7.0 Oct=+3.6
  prcp first5: [5.420, 2.089, 9.031, 7.619, 0.0]  ← matches Julia ✓
  temp first5: [0.228,-0.673, 0.065, 2.729, 1.299] ← Julia reads [-5.986,-6.887,-6.15,-3.485,-4.915]
=#

using Sleipnir, Muninn
using Statistics, Dates, Printf

const RGI = "RGI60-11.01450"
const TSPAN = (2000.0, 2020.0)

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI => rgi_paths[RGI]),
    use_MB = true, use_velocities = false, climate_data_source = :W5E5))
glacier = Sleipnir.initialize_glaciers([RGI], params)[1]

raw = glacier.climate.raw_climate
ref_hgt = glacier.climate.ref_hgt
temp_r = Float64.(vec(raw.temp.data))
prcp_r = Float64.(vec(raw.prcp.data))

println("=== Raw climate RasterStack info ===")
@printf("ref_hgt: %.1f m  n=%d\n", ref_hgt, length(temp_r))
@printf("temp first5:  %s\n", string(round.(temp_r[1:5]; digits = 3)))
@printf("prcp first5:  %s  (Python: [5.420, 2.089, 9.031, 7.619, 0.0])\n",
    string(round.(prcp_r[1:5]; digits = 3)))
@printf("temp mean=%.4f  min=%.2f  max=%.2f\n",
    mean(temp_r), minimum(temp_r), maximum(temp_r))
println("Python: temp mean=+3.9977  min=-17.68  max=+20.76")
println()

# Get the time axis from the raw climate RasterStack
println("=== Time dimension of raw_climate ===")
ti_dim = Sleipnir.dims(raw, Sleipnir.Ti)
ti_vec = collect(ti_dim)
@printf("n times: %d  first=%s  last=%s\n", length(ti_vec), ti_vec[1], ti_vec[end])
println("Python: 7670 days from 1999-01-01 to 2019-12-31")
println()

# Compute monthly means from the raw RasterStack via get_cumulative_climate!
println("=== Monthly means for 2010 via get_cumulative_climate! ===")
println("Month     Julia avg_temp    Python mean")
for m in 1:12
    t_end = Sleipnir.Float(2010.0 + m/12.0)
    step = Sleipnir.Float(1.0/12.0)
    Sleipnir.get_cumulative_climate!(glacier.climate, t_end, step)
    cs = glacier.climate.climate_step
    py = [-6.3, -4.6, -2.1, 2.6, 4.6, 9.8, 13.3, 10.3, 7.0, 3.6, 0.0, -5.6]
    @printf("  2010-%02d:  %8.3f °C    %8.3f °C  (diff=%.3f)\n",
        m, cs.avg_temp, py[m], cs.avg_temp - py[m])
end
println()

# Also compute monthly means directly from the raw data vec (manual)
println("=== Monthly means from raw vec manually ===")
# Get dates from the time axis
dates = Dates.Date.(ti_vec)
for m in 1:12
    idx = findall(d -> year(d)==2010 && month(d)==m, dates)
    if !isempty(idx)
        @printf("  2010-%02d: n=%d  mean=%.3f  min=%.2f  max=%.2f\n",
            m, length(idx), mean(temp_r[idx]),
            minimum(temp_r[idx]), maximum(temp_r[idx]))
    end
end
println()

# Find the index range for August 2010 and print first 5 raw values
println("=== First 5 raw temp values for August 2010 ===")
aug_idx = findall(d -> year(d)==2010 && month(d)==8, dates)
if !isempty(aug_idx)
    @printf("Julia raw values: %s\n",
        string(round.(temp_r[aug_idx[1:min(5, end)]]; digits = 3)))
    @printf("Python values (Aug 2010 first 5 approx): expect ~10°C range\n")
end
