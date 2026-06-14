#=
Decompose default-model MB into accumulation vs melt terms (READ-ONLY probe)
============================================================================
For each glacier, replicate Muninn's daily TI downscaling over the real 2D
ice-cell hypsometry (binned into elevation bands) and report the glacier-wide
area-weighted accumulation term (0.001*prcp_fac*snow) and melt term (DDF*PDD),
under three variants:
  (a) current model   : t_melt = 0,  daily clamped W5E5 lapse rate
  (b) t_melt = -1      : OGGM W5E5 melt threshold
  (c) OGGM lapse       : constant -0.0065 °C/m, t_melt = 0
This shows which single knob flips the default MB from +2 to a realistic value.
=#

using Muninn
using Sleipnir
using Statistics
using Dates
using Printf

const DDF = 0.007    # default TImodel1 DDF (m w.e. /°C/d)
const PRCP_FAC = 1.0      # default
const T_SOLID, T_LIQ = 0.0, 2.0
snow_frac(T) = clamp(1 - (T - T_SOLID) / (T_LIQ - T_SOLID), 0, 1)

glaciers = ["RGI60-11.01450", "RGI60-11.03638"]
tspan = (2000.0, 2020.0)
rgi_paths_all = Sleipnir.get_rgi_paths()

for rgi in glaciers
    rgi_paths = Dict(rgi => rgi_paths_all[rgi])
    params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
        tspan = tspan, multiprocessing = false, workers = 1,
        rgi_paths = rgi_paths, use_MB = true, use_velocities = false,
        climate_data_source = :W5E5))
    glacier = Sleipnir.initialize_glaciers([rgi], params)[1]

    raw = glacier.climate.raw_climate
    ref_hgt = glacier.climate.ref_hgt
    temp = Float64.(collect(raw.temp.data))
    prcp = Float64.(collect(raw.prcp.data))
    grad = clamp.(Float64.(collect(raw.gradient.data)), -0.009, -0.003)
    nyears = tspan[2] - tspan[1]

    # Ice-cell elevations, binned to 25 m for speed (area-weighted by count).
    Svals = glacier.S[.!glacier.mask]
    bins = Dict{Int, Int}()
    for s in Svals
        b = round(Int, s / 25) * 25
        bins[b] = get(bins, b, 0) + 1
    end
    ncells = length(Svals)

    # Accumulate area-weighted snow/PDD (mm and °C·d per year) for 3 variants.
    function decompose(t_melt, gradvec)
        acc_snow = 0.0;
        acc_pdd = 0.0
        for (h, cnt) in bins
            Th = temp .+ gradvec .* (h - ref_hgt)
            snow = sum(prcp .* snow_frac.(Th)) / nyears
            pdd = sum(max.(Th .- t_melt, 0.0)) / nyears
            w = cnt / ncells
            acc_snow += w * snow
            acc_pdd += w * pdd
        end
        acc_term = 0.001 * PRCP_FAC * acc_snow
        melt_term = DDF * acc_pdd
        return acc_term, melt_term
    end

    a_cur, m_cur = decompose(0.0, grad)
    a_tm, m_tm = decompose(-1.0, grad)
    a_la, m_la = decompose(0.0, fill(-0.0065, length(grad)))

    println("="^78)
    println("$rgi   observed geodetic MB = $(round(glacier.geodetic_MB; digits=3)) m/yr")
    println("  ref_hgt=$(round(ref_hgt))  mean(S)=$(round(mean(Svals)))  " *
            "S range $(round(minimum(Svals)))–$(round(maximum(Svals)))  " *
            "lapse(mean)=$(round(mean(grad); digits=5))")
    @printf("  %-22s %10s %10s %10s\n", "variant", "accum", "melt", "net MB")
    @printf("  %-22s %10.3f %10.3f %10.3f\n", "(a) current t_melt=0", a_cur, m_cur,
        a_cur - m_cur)
    @printf("  %-22s %10.3f %10.3f %10.3f\n", "(b) t_melt=-1", a_tm, m_tm, a_tm - m_tm)
    @printf("  %-22s %10.3f %10.3f %10.3f\n", "(c) OGGM lapse -0.0065", a_la, m_la,
        a_la - m_la)
    println()
end
