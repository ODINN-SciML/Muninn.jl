#=
Minimal temperature read test — bypasses glacier initialization
Read the same file at 3 different stages to find where values diverge.
Python truth for raw_climate_(2000.0, 2020.0).nc Aletsch:
  temp first5: [0.228, -0.673, 0.065, 2.729, 1.299]   mean=+3.9977°C
  prcp first5: [5.420, 2.089, 9.031, 7.619, 0.000]
=#

using Sleipnir
using Statistics, Printf

const PATH = joinpath(homedir(),
    ".ODINN/ODINN_prepro/per_glacier/RGI60-11/RGI60-11.01/RGI60-11.01450",
    "raw_climate_(2000.0, 2020.0).nc")

println("File: $PATH")
println("Exists: $(isfile(PATH))")
println()

# ── Stage 1: RasterStack as-is (Float32, no conversion) ──────────────────────
println("Stage 1 — Sleipnir.RasterStack (raw Float32, no conversion):")
rs32 = Sleipnir.RasterStack(PATH)
t32 = vec(rs32.temp.data)
p32 = vec(rs32.prcp.data)
@printf("  temp eltype: %s\n", eltype(t32))
@printf("  temp first5: %s\n", string(round.(Float64.(t32[1:5]); digits = 3)))
@printf("  prcp first5: %s\n", string(round.(Float64.(p32[1:5]); digits = 3)))
@printf("  temp mean=%.4f  min=%.2f  max=%.2f\n",
    mean(Float64.(t32)), minimum(Float64.(t32)), maximum(Float64.(t32)))
println()

# ── Stage 2: After convertRasterStackToFloat64 ────────────────────────────────
println("Stage 2 — after convertRasterStackToFloat64:")
function convertRasterStackToFloat64(rs::Sleipnir.RasterStack)
    layerNames = Sleipnir.Rasters.name(rs)
    return Sleipnir.RasterStack(
        NamedTuple{Tuple(layerNames)}([Float64.(rs[n]) for n in layerNames]),
        metadata = Sleipnir.metadata(rs)
    )
end
rs64 = convertRasterStackToFloat64(rs32)
t64 = vec(rs64.temp.data)
p64 = vec(rs64.prcp.data)
@printf("  temp eltype: %s\n", eltype(t64))
@printf("  temp first5: %s\n", string(round.(t64[1:5]; digits = 3)))
@printf("  prcp first5: %s\n", string(round.(p64[1:5]; digits = 3)))
@printf("  temp mean=%.4f  min=%.2f  max=%.2f\n", mean(t64), minimum(t64), maximum(t64))
println()

println("Python expected:")
println("  temp first5: [0.228, -0.673, 0.065, 2.729, 1.299]  mean=+3.9977")
println("  prcp first5: [5.420, 2.089, 9.031, 7.619, 0.000]")
println()

# ── Stage 3: also check missingval handling ───────────────────────────────────
println("Stage 3 — missingval / fill value info from RasterStack:")
for n in Sleipnir.Rasters.name(rs32)
    r = rs32[n]
    mv = Sleipnir.Rasters.missingval(r)
    @printf("  %s: missingval=%s  first=%s\n", n, mv,
        round(Float64(vec(r.data)[1]); digits = 4))
end

# Count how many values equal or are close to missingval for temp
mv_temp = Sleipnir.Rasters.missingval(rs32[:temp])
n_missing = count(x -> isnan(x) || isinf(x), Float64.(t32))
@printf("  temp: n_inf_or_nan=%d out of %d\n", n_missing, length(t32))
println()

# ── Stage 4: check if the W5E5 original has same issue ───────────────────────
println("Stage 4 — climate_historical_daily_W5E5.nc (original OGGM file):")
w5e5_path = joinpath(homedir(),
    ".ODINN/ODINN_prepro/per_glacier/RGI60-11/RGI60-11.01/RGI60-11.01450",
    "climate_historical_daily_W5E5.nc")
if isfile(w5e5_path)
    rs_w5 = Sleipnir.RasterStack(w5e5_path)
    tw5 = vec(rs_w5.temp.data)
    pw5 = vec(rs_w5.prcp.data)
    @printf("  n records: %d\n", length(tw5))
    @printf("  temp first5: %s\n", string(round.(Float64.(tw5[1:5]); digits = 3)))
    @printf("  prcp first5: %s\n", string(round.(Float64.(pw5[1:5]); digits = 3)))
    @printf("  temp mean=%.4f  min=%.2f  max=%.2f\n",
        mean(Float64.(tw5)), minimum(Float64.(tw5)), maximum(Float64.(tw5)))
    println("  Python (full history, first5): [-8.856, -17.208, -14.391, -12.345, -10.544]")
else
    println("  File not found: $w5e5_path")
end
