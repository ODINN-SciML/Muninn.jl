
export calibrate_ti_model, compute_mean_annual_MB

###############################################
##### GEODETIC MASS BALANCE CALIBRATION  ######
###############################################

"""
    compute_mean_annual_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0),
    ) where {F <: AbstractFloat} -> Sleipnir.Float

Compute the glacier-wide mean annual mass balance (m w.e. yr⁻¹) over the period
`[t_start, t_end]` using a **static** glacier geometry (no ice-flow dynamics).

The surface elevation field `glacier.S` is kept constant throughout the integration.
Climate data is sampled at monthly resolution by default, and the resulting 2D MB
fields are averaged over the ice-covered area (`glacier.mask .== false`) and then
normalised to an annual rate.

As a side-effect this function mutates `glacier.climate` (the cumulative-climate
buffer inside the glacier), which is consistent with the rest of the ODINN stack.

# Arguments

  - `mb_model::TImodel1`: The temperature-index model to evaluate.
  - `glacier::AbstractGlacier`: A glacier object with loaded climate data.
  - `t_start::F`: Start of the integration period (fractional year).
  - `t_end::F`: End of the integration period (fractional year).

# Keyword arguments

  - `step::F`: Integration timestep in fractional years. Default is `1/12` (monthly).

# Returns

  - Glacier-wide mean annual mass balance in m w.e. yr⁻¹.
"""
function compute_mean_annual_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0)) where {F <: AbstractFloat}
    # Accumulate MB over the entire period with static geometry
    total_mb = zeros(Sleipnir.Float, size(glacier.S))

    n_steps = 0
    t = t_start + step
    # Half-step tolerance to guard against floating-point drift at the boundary
    while t <= t_end + step / 2
        get_cumulative_climate!(glacier.climate, Sleipnir.Float(t), Sleipnir.Float(step))
        climate_2D = downscale_2D_climate(
            glacier.climate.climate_step,
            glacier.S,
            glacier.Coords;
            include_topography = false)
        mb_step = compute_MB(mb_model, climate_2D, Sleipnir.Float(step))
        total_mb .+= mb_step
        n_steps += 1
        t += step
    end

    n_steps == 0 && return Sleipnir.Float(0.0)

    # Mask convention: mask == true → outside glacier (no ice)
    glacier_cells = .!glacier.mask
    !any(glacier_cells) && return Sleipnir.Float(0.0)

    # Each monthly compute_MB call returns m w.e. month⁻¹, so the sum over
    # 12 months equals m w.e. yr⁻¹ for one year.  Dividing the total by the
    # number of years gives the mean annual MB.
    n_years = Sleipnir.Float(t_end - t_start)
    cell_values = total_mb[glacier_cells]
    return Sleipnir.Float(sum(cell_values) / length(cell_values)) / n_years
end

"""
    calibrate_ti_model(
        glacier::AbstractGlacier,
        params::Parameters;
        acc_factor::Sleipnir.Float = Sleipnir.Float(1.0 / 1000.0),
      DDF_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
        Sleipnir.Float(params.physical.DDF_min),
        Sleipnir.Float(params.physical.DDF_max)),
      prcp_fac_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
        Sleipnir.Float(params.physical.prcp_fac_min),
        Sleipnir.Float(params.physical.prcp_fac_max)),
        density_ratio::Sleipnir.Float = Sleipnir.Float(1.0),
        calibration_period::Union{Nothing, Tuple{Sleipnir.Float, Sleipnir.Float}} = nothing,
        step::Sleipnir.Float = Sleipnir.Float(1.0 / 12.0),
    ) -> TImodel1

Calibrate `TImodel1` for a single glacier against the geodetic mass balance stored
in `glacier.dhdtData` (e.g. the 2000-2020 observations from Hugonnet et al. 2021).

The calibration follows a **2-step cascade** analogous to the OGGM v1.6 approach,
using Brent's method for root-finding at each step:

1. **DDF step**: With `prcp_fac = 1.0` fixed, find the degree-day factor that
   matches the geodetic MB target.  If the target can be bracketed within
   `DDF_bounds`, this step alone produces the calibrated model.

2. **`prcp_fac` step** (fallback): If the DDF search cannot bracket the target
   (e.g. because the climate input severely under/overestimates precipitation),
   the DDF is fixed at its boundary value and `prcp_fac` is varied within
   `prcp_fac_bounds` instead.  A warning is emitted when this fallback is used.

A static glacier geometry is assumed throughout (no ice-flow dynamics).

# Arguments

  - `glacier::AbstractGlacier`: Glacier whose `dhdtData` field contains the
    observed geodetic mass balance (see [`DhdtData`](@ref)).
  - `params::Parameters`: Simulation parameters (used only for type consistency).

# Keyword arguments

  - `acc_factor`: Unit-conversion accumulation factor (m w.e. mm⁻¹), kept fixed.
    Default: `1.0 × 10⁻³`.
  - `DDF_bounds`: Search interval `(DDF_min, DDF_max)` in m w.e. °C⁻¹ d⁻¹.
    Default: `(0.5 × 10⁻³, 20.0 × 10⁻³)`.
  - `prcp_fac_bounds`: Search interval for the precipitation correction factor
    (dimensionless). Default: `(0.1, 10.0)`.
  - `density_ratio`: Conversion factor applied to `glacier.dhdtData.dhdt`.  Use
    `params.physical.ρ / params.physical.ρ_w ≈ 0.9` when the data are in m ice yr⁻¹, or `1.0`
    (default) for m w.e. yr⁻¹ (Hugonnet et al. 2021).
  - `calibration_period`: Time window `(t_start, t_end)` in fractional years.
    Defaults to `glacier.dhdtData.t`.
  - `step`: Integration timestep in fractional years. Default: `1/12` (monthly).

# Returns

  - A `TImodel1{Sleipnir.Float}` with calibrated `DDF` and `prcp_fac`.

# Reference

Hugonnet, R. et al. (2021). Accelerated global glacier mass loss in the early
twenty-first century. *Nature*, 592, 726–731.
https://doi.org/10.1038/s41586-021-03436-z
"""
function calibrate_ti_model(
        glacier::AbstractGlacier,
        params::Parameters;
        acc_factor::Sleipnir.Float = Sleipnir.Float(1.0 / 1000.0),
        DDF_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
            Sleipnir.Float(params.physical.DDF_min),
            Sleipnir.Float(params.physical.DDF_max)),
        prcp_fac_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
            Sleipnir.Float(params.physical.prcp_fac_min),
            Sleipnir.Float(params.physical.prcp_fac_max)),
        density_ratio::Sleipnir.Float = Sleipnir.Float(1.0),
        calibration_period::Union{Nothing, Tuple{Sleipnir.Float, Sleipnir.Float}} = nothing,
        step::Sleipnir.Float = Sleipnir.Float(1.0 / 12.0))
    if isnothing(glacier.dhdtData)
        throw(ArgumentError(
            "glacier.dhdtData is nothing. Geodetic mass-balance observations " *
            "are required for TI model calibration. " *
            "Please populate glacier.dhdtData (e.g. from Hugonnet et al. 2021) " *
            "before calling calibrate_ti_model."))
    end

    # Determine calibration period
    t_start, t_end = if isnothing(calibration_period)
        (Sleipnir.Float(glacier.dhdtData.t[1]),
            Sleipnir.Float(glacier.dhdtData.t[2]))
    else
        calibration_period
    end

    # Target glacier-wide mean annual MB in m w.e. yr⁻¹
    mb_target = Sleipnir.Float(glacier.dhdtData.dhdt) * density_ratio

    # ── Step 1: calibrate DDF with prcp_fac = 1.0 ──────────────────────────
    prcp_fac_fixed = Sleipnir.Float(1.0)

    function residual_ddf(DDF::Sleipnir.Float)
        model_trial = TImodel1{Sleipnir.Float}(DDF, acc_factor, prcp_fac_fixed)
        return compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) - mb_target
    end

    DDF_min, DDF_max = DDF_bounds
    r_ddf_min = residual_ddf(DDF_min)
    r_ddf_max = residual_ddf(DDF_max)

    if r_ddf_min * r_ddf_max <= 0
        # Root is bracketed → Brent solve
        DDF_cal = _brent(residual_ddf, DDF_min, DDF_max, r_ddf_min, r_ddf_max)
        return TImodel1{Sleipnir.Float}(DDF_cal, acc_factor, prcp_fac_fixed)
    end

    # ── Step 2: DDF at boundary, calibrate prcp_fac ────────────────────────
    # Fix DDF at whichever boundary is closer to the target
    DDF_fixed = abs(r_ddf_min) <= abs(r_ddf_max) ? DDF_min : DDF_max

    @warn "calibrate_ti_model: geodetic MB target ($(round(mb_target; digits=4)) m w.e. yr⁻¹) " *
          "for glacier $(glacier.rgi_id) could not be bracketed by DDF alone within " *
          "DDF_bounds = $(DDF_bounds). " *
          "Falling back to prcp_fac calibration with DDF fixed at $(DDF_fixed)."

    function residual_prcp(prcp_fac::Sleipnir.Float)
        model_trial = TImodel1{Sleipnir.Float}(DDF_fixed, acc_factor, prcp_fac)
        return compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) - mb_target
    end

    pf_min, pf_max = prcp_fac_bounds
    r_pf_min = residual_prcp(pf_min)
    r_pf_max = residual_prcp(pf_max)

    if r_pf_min * r_pf_max <= 0
        prcp_fac_cal = _brent(residual_prcp, pf_min, pf_max, r_pf_min, r_pf_max)
        return TImodel1{Sleipnir.Float}(DDF_fixed, acc_factor, prcp_fac_cal)
    end

    # Both steps exhausted — return the boundary pair with the smallest residual
    @warn "calibrate_ti_model: prcp_fac calibration also failed to bracket the target " *
          "for glacier $(glacier.rgi_id). Returning the best boundary pair."
    prcp_fac_cal = abs(r_pf_min) <= abs(r_pf_max) ? pf_min : pf_max
    return TImodel1{Sleipnir.Float}(DDF_fixed, acc_factor, prcp_fac_cal)
end

"""
    calibrate_ti_model(
        glaciers::Vector{<:AbstractGlacier},
        params::Parameters;
        kwargs...,
    ) -> Vector{TImodel1}

Calibrate `TImodel1` independently for each glacier in `glaciers`.

All keyword arguments are forwarded to the single-glacier method; see
[`calibrate_ti_model(glacier, params; ...)`](@ref) for details.

Glaciers whose `dhdtData` is `nothing` are skipped and a warning is emitted.
Returns a `Vector{TImodel1}` of the same length as `glaciers`.  Entries
corresponding to skipped glaciers contain the default (uncalibrated) model.
"""
function calibrate_ti_model(
        glaciers::Vector{<:AbstractGlacier},
        params::Parameters;
        kwargs...)
    map(glaciers) do glacier
        if isnothing(glacier.dhdtData)
            @warn "calibrate_ti_model: skipping glacier $(glacier.rgi_id) " *
                  "because dhdtData is nothing."
            return TImodel1(params)
        end
        return calibrate_ti_model(glacier, params; kwargs...)
    end
end

# ---- private root-finder -----------------------------------------------

# Thin wrapper around Roots.find_zero so call-sites stay unchanged.
# Uses Brent's method (same algorithm as scipy.optimize.brentq in OGGM).
function _brent(f, a::F, b::F, fa::F, fb::F;
        tol::F = F(1e-9), max_iter::Int = 100) where {F <: AbstractFloat}
    return F(find_zero(f, (a, b), Brent();
        atol = tol, rtol = zero(F), maxiters = max_iter))
end
