
export calibrate_MB_model, calibrate_ti_model, compute_mean_annual_MB,
       compute_cumulative_MB

###############################################
##### GEODETIC MASS BALANCE CALIBRATION  ######
###############################################

"""
    compute_cumulative_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0),
    ) where {F <: AbstractFloat} -> Matrix{Sleipnir.Float}

Accumulate the gridded mass balance (m w.e.) over `[t_start, t_end]` using a
static glacier geometry.  Returns the full 2D cumulative MB field, suitable for
spatial visualisation via [`plot_cumulative_mb`](@ref).

See [`compute_mean_annual_MB`](@ref) for the glacier-wide scalar summary.
"""
function compute_cumulative_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0)) where {F <: AbstractFloat}
    total_mb = zeros(Sleipnir.Float, size(glacier.S))
    t = t_start + step
    while t <= t_end + step / 2
        get_cumulative_climate!(glacier.climate, Sleipnir.Float(t), Sleipnir.Float(step))
        climate_2D = downscale_2D_climate(
            glacier.climate.climate_step,
            glacier.climate.climate_raw_step,
            glacier.S,
            glacier.Coords;
            include_topography = false,
            temp_bias = mb_model.temp_bias)
        total_mb .+= compute_MB(mb_model, climate_2D, Sleipnir.Float(step))
        t += step
    end
    return total_mb
end

"""
    compute_mean_annual_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0),
    ) where {F <: AbstractFloat} -> Sleipnir.Float

Glacier-wide mean annual mass balance (m w.e. yr⁻¹) over `[t_start, t_end]`
with static geometry.  Spatially averages the output of [`compute_cumulative_MB`](@ref)
over the ice-covered area and normalises to an annual rate.
"""
function compute_mean_annual_MB(
        mb_model::TImodel1,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0)) where {F <: AbstractFloat}
    total_mb = compute_cumulative_MB(mb_model, glacier, t_start, t_end; step)
    # Empty calibration window: compute_cumulative_MB returns a zero matrix,
    # detected below via the glacier mask / iszero check.
    glacier_cells = .!glacier.mask
    !any(glacier_cells) && return Sleipnir.Float(0.0)
    n_years = Sleipnir.Float(t_end - t_start)
    iszero(total_mb) && return Sleipnir.Float(0.0)
    cell_values = total_mb[glacier_cells]
    return Sleipnir.Float(sum(cell_values) / length(cell_values)) / n_years
end

"""
    calibrate_ti_model(
        glacier::AbstractGlacier,
        params::Parameters;
        DDF_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
            Sleipnir.Float(params.physical.DDF_min),
            Sleipnir.Float(params.physical.DDF_max)),
        prcp_fac_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
            Sleipnir.Float(params.physical.prcp_fac_min),
            Sleipnir.Float(params.physical.prcp_fac_max)),
        temp_bias_bounds::Tuple{Sleipnir.Float, Sleipnir.Float} = (
            Sleipnir.Float(params.physical.temp_bias_min),
            Sleipnir.Float(params.physical.temp_bias_max)),
        density_ratio::Sleipnir.Float = Sleipnir.Float(1.0),
        calibration_period::Union{Nothing, Tuple{Sleipnir.Float, Sleipnir.Float}} = nothing,
        prcp_fac::Union{Symbol, Real} = :from_winter_prcp,
        step::Sleipnir.Float = Sleipnir.Float(1.0 / 12.0),
    ) -> TImodel1

Calibrate `TImodel1` for a single glacier against the geodetic mass balance stored
in `glacier.dhdtData` (e.g. the 2000-2020 observations from Hugonnet et al. 2021).

The calibration follows a **3-step cascade** analogous to the OGGM v1.6 approach,
using Brent's method for root-finding at each step:

 1. **DDF step**: With `prcp_fac` fixed (glacier-specific from winter precipitation
    by default, see the `prcp_fac` keyword) and `temp_bias = 0.0`, find the
    degree-day factor that matches the geodetic MB target.  If the target can be
    bracketed within `DDF_bounds`, this step alone produces the calibrated model.

 2. **`prcp_fac` step** (fallback): If the DDF search cannot bracket the target,
    DDF is fixed at its boundary value and `prcp_fac` is varied within
    `prcp_fac_bounds`.  A warning is emitted when this fallback is used.

 3. **`temp_bias` step** (fallback): If both previous steps fail, DDF and
    `prcp_fac` are fixed at their best boundary values and a uniform temperature
    bias (°C) is varied within `temp_bias_bounds`.  This handles glaciers where
    the climate forcing is systematically biased for the glacier's hypsometry
    (e.g. high-elevation accumulation overestimation in coarse reanalysis data).

A static glacier geometry is assumed throughout (no ice-flow dynamics).

# Arguments

  - `glacier::AbstractGlacier`: Glacier whose `dhdtData` field contains the
    observed geodetic mass balance (see [`DhdtData`](@ref)).
  - `params::Parameters`: Simulation parameters (provides physical bounds).

# Keyword arguments

  - `DDF_bounds`: Search interval `(DDF_min, DDF_max)` in m w.e. °C⁻¹ d⁻¹.
    Default: from `params.physical`.
  - `prcp_fac_bounds`: Search interval for the precipitation correction factor
    (dimensionless). Default: from `params.physical`.
  - `temp_bias_bounds`: Search interval for the temperature bias (°C).
    Default: from `params.physical`.
  - `density_ratio`: Conversion factor applied to `glacier.dhdtData.dhdt`.  Use
    `params.physical.ρ / params.physical.ρ_w ≈ 0.9` when the data are in m ice yr⁻¹, or `1.0`
    (default) for m w.e. yr⁻¹ (Hugonnet et al. 2021).
  - `calibration_period`: Time window `(t_start, t_end)` in fractional years.
    Defaults to `glacier.dhdtData.t`.
  - `prcp_fac`: Precipitation factor used in the DDF step. `:from_winter_prcp`
    (default) derives a glacier-specific factor from winter precipitation via
    [`Sleipnir.get_winter_prcp_factor`](@ref); a `Real` value fixes it (e.g. `2.5`
    for OGGM's global W5E5 default). Only used in step 1; the fallback steps still
    search `prcp_fac_bounds`.
  - `step`: Integration timestep in fractional years. Default: `1/12` (monthly).

# Returns

  - A `TImodel1{Sleipnir.Float}` with calibrated `DDF`, `prcp_fac`, and `temp_bias`.

This is the per-glacier building block used by [`calibrate_MB_model`](@ref) to
build a per-glacier vector of calibrated models.

# References

Hugonnet, R. et al. (2021). Accelerated global glacier mass loss in the early
twenty-first century. *Nature*, 592, 726–731.
https://doi.org/10.1038/s41586-021-03436-z

Maussion, F., Butenko, A., Eis, J., Fourteau, K., Jarosch, A. H., Landmann, J.,
Oesterle, F., Recinos, B., USias, S., Valsecchi, L., Marzeion, B., and Cogley,
J. G. (2019). The Open Global Glacier Model (OGGM) v1.1. *Geoscientific Model
Development*, 12, 909–941.
https://doi.org/10.5194/gmd-12-909-2019
"""
function calibrate_ti_model(
        glacier::G,
        params::Parameters;
        DDF_bounds::Tuple{<:AbstractFloat, <:AbstractFloat} = (
            params.physical.DDF_min, params.physical.DDF_max),
        prcp_fac_bounds::Tuple{<:AbstractFloat, <:AbstractFloat} = (
            params.physical.prcp_fac_min, params.physical.prcp_fac_max),
        temp_bias_bounds::Tuple{<:AbstractFloat, <:AbstractFloat} = (
            params.physical.temp_bias_min, params.physical.temp_bias_max),
        density_ratio::Real = 1.0,
        calibration_period::Union{Nothing, Tuple{<:AbstractFloat, <:AbstractFloat}} = nothing,
        prcp_fac::Union{Symbol, Real} = :from_winter_prcp,
        step::Real = 1.0 / 12.0) where {G <: AbstractGlacier}
    if isnothing(glacier.dhdtData) || !isfinite(glacier.dhdtData.dhdt)
        throw(ArgumentError(
            "glacier.dhdtData.dhdt is not available for $(glacier.rgi_id). " *
            "Geodetic mass-balance observations are required for TI model calibration. " *
            "Please populate glacier.dhdtData (e.g. from Hugonnet et al. 2021) " *
            "before calling calibrate_ti_model."))
    end
    mb_observation = glacier.dhdtData.dhdt

    # Determine calibration period
    t_start,
    t_end = if isnothing(calibration_period)
        (Sleipnir.Float(glacier.dhdtData.t[1]),
            Sleipnir.Float(glacier.dhdtData.t[2]))
    else
        calibration_period
    end

    # Target glacier-wide mean annual MB in m w.e. yr⁻¹
    mb_target = Sleipnir.Float(mb_observation) * density_ratio

    tb_zero = Sleipnir.Float(0.0)

    # ── Step 1: calibrate DDF with prcp_fac fixed, temp_bias = 0.0 ────────
    # prcp_fac is either glacier-specific (from winter precipitation, OGGM-style)
    # or a user-provided constant (e.g. 2.5, OGGM's global W5E5 default).
    prcp_fac_fixed = if prcp_fac === :from_winter_prcp
        Sleipnir.get_winter_prcp_factor(glacier, params; prcp_fac_bounds)
    else
        Sleipnir.Float(prcp_fac)
    end

    function residual_ddf(DDF::Sleipnir.Float)
        model_trial = TImodel1{Sleipnir.Float}(DDF, prcp_fac_fixed, tb_zero)
        return compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) -
               mb_target
    end

    DDF_min, DDF_max = DDF_bounds
    r_ddf_min = residual_ddf(DDF_min)
    r_ddf_max = residual_ddf(DDF_max)

    if r_ddf_min * r_ddf_max <= 0
        DDF_cal = _brent(residual_ddf, DDF_min, DDF_max, r_ddf_min, r_ddf_max)
        return TImodel1{Sleipnir.Float}(DDF_cal, prcp_fac_fixed, tb_zero)
    end

    # ── Step 2: DDF at boundary, calibrate prcp_fac ────────────────────────
    DDF_fixed = abs(r_ddf_min) <= abs(r_ddf_max) ? DDF_min : DDF_max

    @warn "calibrate_ti_model: geodetic MB target ($(round(mb_target; digits=4)) m w.e. yr⁻¹) " *
          "for glacier $(glacier.rgi_id) could not be bracketed by DDF alone within " *
          "DDF_bounds = $(DDF_bounds). " *
          "Falling back to prcp_fac calibration with DDF fixed at $(DDF_fixed)."

    function residual_prcp(prcp_fac::Sleipnir.Float)
        model_trial = TImodel1{Sleipnir.Float}(DDF_fixed, prcp_fac, tb_zero)
        return compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) -
               mb_target
    end

    pf_min, pf_max = prcp_fac_bounds
    r_pf_min = residual_prcp(pf_min)
    r_pf_max = residual_prcp(pf_max)

    if r_pf_min * r_pf_max <= 0
        prcp_fac_cal = _brent(residual_prcp, pf_min, pf_max, r_pf_min, r_pf_max)
        return TImodel1{Sleipnir.Float}(DDF_fixed, prcp_fac_cal, tb_zero)
    end

    # ── Step 3: DDF + prcp_fac at boundary, calibrate temp_bias ───────────
    prcp_fac_fixed2 = abs(r_pf_min) <= abs(r_pf_max) ? pf_min : pf_max

    @warn "calibrate_ti_model: prcp_fac calibration also failed to bracket the target " *
          "for glacier $(glacier.rgi_id). " *
          "Falling back to temp_bias calibration with DDF=$(DDF_fixed), prcp_fac=$(prcp_fac_fixed2)."

    function residual_tb(temp_bias::Sleipnir.Float)
        model_trial = TImodel1{Sleipnir.Float}(
            DDF_fixed, prcp_fac_fixed2, temp_bias)
        return compute_mean_annual_MB(model_trial, glacier, t_start, t_end; step) -
               mb_target
    end

    tb_min, tb_max = temp_bias_bounds
    r_tb_min = residual_tb(tb_min)
    r_tb_max = residual_tb(tb_max)

    if r_tb_min * r_tb_max <= 0
        tb_cal = _brent(residual_tb, tb_min, tb_max, r_tb_min, r_tb_max)
        return TImodel1{Sleipnir.Float}(DDF_fixed, prcp_fac_fixed2, tb_cal)
    end

    # All steps exhausted — return the boundary triple with the smallest residual
    @warn "calibrate_ti_model: temp_bias calibration also failed to bracket the target " *
          "for glacier $(glacier.rgi_id). Returning the best boundary triple."
    tb_cal = abs(r_tb_min) <= abs(r_tb_max) ? tb_min : tb_max
    return TImodel1{Sleipnir.Float}(DDF_fixed, prcp_fac_fixed2, tb_cal)
end

"""
    calibrate_MB_model(
        model::Sleipnir.Model,
        glaciers::Vector{<:AbstractGlacier},
        params::Parameters,
    ) -> Sleipnir.Model

High-level entry point that calibrates the mass balance model of `model`
per glacier against geodetic observations, returning a (possibly new) `model`.

Calibration cannot happen in place: for `TImodel1` it replaces a single model by one
model per glacier, which changes the type of the `mass_balance` field. Callers must use
the return value; discarding it silently keeps the uncalibrated model.

The behaviour dispatches on the mass balance model type:

  - `TImodel1`: returns a new `Model` whose `mass_balance` is a per-glacier vector of
    `TImodel1`s, each fitted against its glacier's geodetic mass balance with
    [`calibrate_ti_model`](@ref).  Glaciers without `dhdtData` keep the original
    (uncalibrated) model and a warning is emitted.  If no glacier carries geodetic
    data, the model is returned unchanged.
  - any other [`MBmodel`](@ref): no-op — no calibration routine is defined, so the
    model is returned untouched.  Add a `_calibrate_MB_model(model, ::YourType, …)`
    method to support a new model type.

An already-vectorized (per-glacier) mass balance model is left unchanged.  Any
keyword arguments are forwarded to the type-specific calibrator (e.g.
[`calibrate_ti_model`](@ref) for `TImodel1`).  This is the function the
`Prediction` and `Inversion` constructors call when
`params.simulation.calibrate_MB` is `true`.
"""
function calibrate_MB_model(
        model::Sleipnir.Model,
        glaciers::Vector{<:AbstractGlacier},
        params::Parameters;
        kwargs...)
    # Nothing to do for empty or already per-glacier mass balance models.
    isnothing(model.mass_balance) && return model
    model.mass_balance isa AbstractVector && return model
    return _calibrate_MB_model(model, model.mass_balance, glaciers, params; kwargs...)
end

# Default: mass balance models without a calibration routine are left untouched.
function _calibrate_MB_model(
        model, ::MBmodel, ::Vector{<:AbstractGlacier}, ::Parameters; kwargs...)
    return model
end

# TImodel1: fit one model per glacier against its geodetic mass balance.
# The outer function is parametric in G so the pmap closure captures G as a
# concrete type parameter. The `glacier::G` typeassert inside the closure
# narrows pmap's Channel{Any} element to G, keeping calibrate_ti_model
# dispatch and _brent fully type-stable for JET.
function _calibrate_MB_model(
        model, template::TImodel1,
        glaciers::Vector{G}, params::Parameters;
        kwargs...) where {G <: AbstractGlacier}
    # Without any geodetic observations there is nothing to calibrate against.
    any(g -> !isnothing(g.dhdtData), glaciers) || return model
    results = pmap(glaciers) do glacier_any
        glacier = glacier_any::G
        if isnothing(glacier.dhdtData)
            @warn "calibrate_MB_model: skipping glacier $(glacier.rgi_id) " *
                  "because dhdtData is nothing."
            return template
        end
        return calibrate_ti_model(glacier, params; kwargs...)
    end
    calibrated = _narrow_models(typeof(template), results)
    return Sleipnir.Model(model.iceflow, calibrated, model.trainable_components)
end

# Narrow the pmap result (Vector{Any} via Channel{Any}) to the concrete model
# type via per-element typeasserts so the rebuilt Model stays type-stable.
function _narrow_models(::Type{T}, results::AbstractVector) where {T}
    out = Vector{T}(undef, length(results))
    @inbounds for i in eachindex(results)
        out[i] = results[i]::T
    end
    return out
end

# ---- private root-finder -----------------------------------------------

# Thin wrapper around Roots.find_zero so call-sites stay unchanged.
# Uses Brent's method (same algorithm as scipy.optimize.brentq in OGGM).
function _brent(f, a::F, b::F, fa::F, fb::F;
        tol::F = F(1e-9), max_iter::Int = 100) where {F <: AbstractFloat}
    return F(find_zero(f, (a, b), Brent();
        atol = tol, rtol = zero(F), maxiters = max_iter))
end
