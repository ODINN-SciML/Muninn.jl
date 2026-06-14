using Muninn
using Sleipnir
using Statistics
using Dates

rgi = "RGI60-11.03638"  # Argentière
tspan = (2000.0, 2020.0)
rgi_paths = Sleipnir.get_rgi_paths()
rgi_paths = Dict(rgi => rgi_paths[rgi])
params = Sleipnir.Parameters(
    simulation = Sleipnir.SimulationParameters(
    tspan = tspan, multiprocessing = false, workers = 1,
    rgi_paths = rgi_paths, use_MB = true, use_velocities = false))
glacier = Sleipnir.initialize_glaciers([rgi], params)[1]

clim = glacier.climate
raw = clim.raw_climate
ref_hgt = clim.ref_hgt
println("ref_hgt = ", ref_hgt)

# Daily raw series over the whole record
temp = Float64.(collect(raw.temp.data))       # daily temperature at ref_hgt (°C)
prcp = Float64.(collect(raw.prcp.data))        # daily precipitation (mm)
grad = Float64.(collect(raw.gradient.data))    # daily lapse rate (°C/m)
grad = clamp.(grad, -0.009, -0.003)
dates = collect(Sleipnir.dims(raw, Sleipnir.Ti))
ym = [(year(d), month(d)) for d in dates]
nyears = (tspan[2] - tspan[1])
println("n daily records = ", length(temp), "  over ",
    round(length(temp)/365.25; digits = 1), " yr")

t_solid, t_liq = 0.0, 2.0
snow_frac(T) = clamp(1 - (T - t_solid)/(t_liq - t_solid), 0, 1)

for h in (1600.0, 2200.0, 2800.0, 3400.0)
    Th = temp .+ grad .* (h - ref_hgt)                 # daily temp at elevation h

    # --- accumulation: solid precip ---
    snow_daily = sum(prcp .* snow_frac.(Th)) / nyears  # OGGM daily linear ramp
    # monthly hard threshold (our current model): all month's prcp is snow if monthly mean <= 0
    snow_monthly = 0.0
    for k in unique(ym)
        idx = findall(==(k), ym)
        if mean(Th[idx]) <= 0.0
            snow_monthly += sum(prcp[idx])
        end
    end
    snow_monthly /= nyears

    # --- melt: degree days ---
    dd_tmelt0 = sum(max.(Th, 0.0)) / nyears           # ours (t_melt = 0)
    dd_tmelt_1 = sum(max.(Th .+ 1.0, 0.0)) / nyears     # OGGM (t_melt = -1)

    println("\nelevation $(Int(h)) m:")
    println("  annual solid precip:  daily(OGGM)=", round(snow_daily; digits = 0),
        " mm   monthly-hard(ours)=", round(snow_monthly; digits = 0), " mm")
    println("  annual degree-days:   t_melt=0(ours)=", round(dd_tmelt0; digits = 0),
        "   t_melt=-1(OGGM)=", round(dd_tmelt_1; digits = 0))
end
