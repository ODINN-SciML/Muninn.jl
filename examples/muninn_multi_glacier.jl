#=
Multi-glacier Muninn calibration — write per-cell MB stats to JSON
====================================================================
Calibrates TImodel1 for 7 WGMS glaciers against Hugonnet 2000-2020,
computes per 2D-cell annual MB statistics (m w.e. yr⁻¹), and writes
muninn_stats.json for the Python compare_argentiere_oggm.py script.
=#

using Muninn, Sleipnir, Statistics, DelimitedFiles, Printf

const TSPAN = (2000.0, 2020.0)
const NYEARS = Float64(TSPAN[2] - TSPAN[1])

const GLACIERS = [
    "RGI60-11.03638",  # Argentière
    "RGI60-11.01450",  # Aletschgletscher
    "RGI60-11.01238",  # Rhonegletscher
    "RGI60-11.00897",  # Hintereisferner
    "RGI60-08.00213",  # Storglaciaeren
    "RGI60-01.00570",  # Gulkana
    "RGI60-02.05098"  # Peyto
]

const NAMES = Dict(
    "RGI60-11.03638" => "Argentière",
    "RGI60-11.01450" => "Aletschgletscher",
    "RGI60-11.01238" => "Rhonegletscher",
    "RGI60-11.00897" => "Hintereisferner",
    "RGI60-08.00213" => "Storglaciaeren",
    "RGI60-01.00570" => "Gulkana",
    "RGI60-02.05098" => "Peyto"
)

rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(r => rgi_paths[r] for r in GLACIERS),
    use_MB = true, use_velocities = false, climate_data_source = :W5E5))

println("Initializing $(length(GLACIERS)) glaciers...")
glaciers = Sleipnir.initialize_glaciers(GLACIERS, params)
println("Done.\n")

rows = []  # collect rows for CSV

for glacier in glaciers
    rgi = glacier.rgi_id
    name = NAMES[rgi]
    println("=" ^ 60)
    println("$name  ($rgi)")

    # Calibrate
    cal = calibrate_ti_model(glacier, params)

    t0 = Sleipnir.Float(TSPAN[1])
    t1 = Sleipnir.Float(TSPAN[2])

    # Per-cell cumulative MB → annual (m w.e. yr⁻¹)
    total_mb = Muninn.compute_cumulative_MB(cal, glacier, t0, t1)
    ann_mb = total_mb ./ NYEARS        # m w.e. yr⁻¹ per cell

    # Ice cells only
    ice_mask = .!glacier.mask
    cells = ann_mb[ice_mask]
    n_cells = count(ice_mask)

    mf = Float64(cal.DDF * 1000)
    pf = Float64(cal.prcp_fac)
    tb = Float64(get_temp_bias(cal))
    mn = Float64(mean(cells))
    med = Float64(median(cells))
    sd = Float64(std(cells))
    lo = Float64(minimum(cells))
    hi = Float64(maximum(cells))

    area_km2 = count(ice_mask) * glacier.Δx * glacier.Δy / 1e6
    push!(rows, [rgi, name, area_km2, glacier.geodetic_MB, n_cells,
        mf, pf, tb, mn, med, sd, lo, hi])

    @printf("  melt_f=%.4f mm/d/K  pf=%.4f  tb=%.4f  n_cells=%d\n", mf, pf, tb, n_cells)
    @printf("  MB stats: mean=%.3f  med=%.3f  std=%.3f  min=%.3f  max=%.3f\n\n",
        mn, med, sd, lo, hi)
end

# Write CSV for Python to read
out_path = joinpath(
    @__DIR__, "..", "..", "..", "..", "Python", "Gungnir", "muninn_stats.csv")
header = ["rgi_id", "name", "area_km2", "hugonnet", "n_cells",
    "melt_f", "prcp_fac", "temp_bias", "mean", "median", "std", "min", "max"]
open(out_path, "w") do io
    println(io, join(header, ","))
    for row in rows
        println(io, join(row, ","))
    end
end
println("Wrote: $out_path")
