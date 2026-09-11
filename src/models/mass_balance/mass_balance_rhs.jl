export MBcache, ElevationLUT, MB_rate!, MB_rate_∂H!, MB_rate_∂H_maxabs,
       build_elevation_lut, mb_S_dependence, smoothstep, smoothstep_∂, mb_cache_active

import Sleipnir: init_mb_cache, mb_cache_type

# Ice thickness scales (m) that gate where the mass balance rate acts.
#
# Ablation ramps from zero over `H_ABL_DEFAULT`, so a cell can't melt through the bed.
#
# Accumulation switches on centred on `H_ACC_DEFAULT`, over a width `H_ACC_W_DEFAULT` —
# not ramped from zero. A marginal accumulation cell must reach the full rate before it can
# thicken enough to spill ice over a divide into the neighbouring catchment; ramping from
# zero would feed such cells at half rate long before that.
const H_ACC_DEFAULT = Sleipnir.Float(10.0)
const H_ACC_W_DEFAULT = Sleipnir.Float(1.0)
const H_ABL_DEFAULT = Sleipnir.Float(1.0)

# Elevation resolution (m) of the mass balance lookup table, and the padding (m) added
# either side of the elevation range the glacier can reach during a run.
const LUT_dΔS_DEFAULT = Sleipnir.Float(1.0)
const LUT_PAD_DEFAULT = Sleipnir.Float(1000.0)

"""
    mb_S_dependence(mb_model::MBmodel)

How the mass balance rate of `mb_model` depends on the ice surface `S`, which decides how
it can be evaluated inside the ice flow RHS.

  - `:elevation_only`: the rate at a cell depends on `S` only through the scalar offset
    `ΔS = S - ref_hgt`. The whole model then collapses to a one-dimensional function of
    `ΔS` per time window, which is what makes a lookup table possible.
  - `:general`: the rate depends on `S` in a way that does not reduce to `ΔS` (a neural
    network reading several fields, for instance) and must be evaluated directly.

`:general` is the conservative default; models opt in to `:elevation_only`.
"""
mb_S_dependence(::MBmodel) = :general
mb_S_dependence(::TImodel1) = :elevation_only

"""
    smoothstep(x)

Cubic Hermite ramp `x²(3 - 2x)`, clamped to `[0, 1]`.

Fades the mass balance rate in and out with ice thickness. A hard mask (as
`apply_MB_mask!` applies) is a step discontinuity in `H`, which as a source term is
unresolvable by an adaptive controller and undifferentiable for an adjoint. This ramp is C¹,
with derivative bounded by `1.5` over its width.
"""
@inline function smoothstep(x::R) where {R <: Real}
    x <= zero(R) && return zero(R)
    x >= one(R) && return one(R)
    return x * x * (R(3) - R(2) * x)
end

"""
    smoothstep_∂(x)

Derivative of [`smoothstep`](@ref), zero outside `[0, 1]`.
"""
@inline function smoothstep_∂(x::R) where {R <: Real}
    (x <= zero(R) || x >= one(R)) && return zero(R)
    return R(6) * x * (one(R) - x)
end

"""
    ElevationLUT{F <: AbstractFloat}

Lookup table of positive degree days and solid precipitation as a function of the elevation
offset `ΔS = S - ref_hgt`, one column per mass balance window.

Both quantities are *exactly* piecewise linear in `ΔS`:

```
PDD(ΔS)  = Σ_d max(0, T_d + g_d·ΔS)
snow(ΔS) = Σ_d p_d·clamp((2 - T_d - g_d·ΔS)/2, 0, 1)
```

so linear interpolation on a uniform table is exact except within one cell of a breakpoint.
Monthly mass balance error from the table is 7.1e-5 m at `dΔS = 1` m; the lookup costs about
an eighth of one `SIA2D!` call, versus ten times that for the daily loop it replaces.

# Fields

  - `ΔS_min::F`, `ΔS_max::F`: Elevation offset range covered, in m.
  - `dΔS::F`, `inv_dΔS::F`: Table resolution in m, and its reciprocal.
  - `PDD::Matrix{F}`: Positive degree days, `(n_elevations, n_windows)`.
  - `snow::Matrix{F}`: Solid precipitation in mm, same shape.
"""
struct ElevationLUT{F <: AbstractFloat}
    ΔS_min::F
    ΔS_max::F
    dΔS::F
    inv_dΔS::F
    PDD::Matrix{F}
    snow::Matrix{F}
end

Base.isempty(lut::ElevationLUT) = isempty(lut.PDD)

"""
    pdd_snow_exact(window::ClimateWindow, ΔS, temp_bias)

Exact daily sums of positive degree days and solid precipitation at elevation offset `ΔS`.

This is the reference the lookup table approximates, and the routine that fills it. It is
the same daily temperature-index computation `downscale_2D_climate!` performs, evaluated at
a single elevation instead of over a grid.
"""
@inline function pdd_snow_exact(
        window::ClimateWindow, ΔS::R, temp_bias::Real) where {R <: Real}
    pdd = zero(R)
    snow = zero(R)
    @inbounds for d in eachindex(window.temp)
        T = window.temp[d] + temp_bias
        g = window.gradient[d]
        pdd += max(zero(R), T + g * ΔS)
        snow += window.prcp[d] * clamp((2 - T - g * ΔS) / 2, zero(R), one(R))
    end
    return pdd, snow
end

"""
    build_elevation_lut(windows, ΔS_min, ΔS_max, temp_bias; dΔS = LUT_dΔS_DEFAULT)

Build an [`ElevationLUT`](@ref) covering `[ΔS_min, ΔS_max]` for the given climate windows.

`temp_bias` is baked into the table, so a table is only valid for the mass balance model it
was built from. `DDF` and `prcp_fac` are deliberately *not* baked in: they scale the melt
and accumulation terms linearly and are applied at lookup time, which keeps the table valid
if they change.
"""
function build_elevation_lut(
        windows::AbstractVector{<:ClimateWindow},
        ΔS_min::Real, ΔS_max::Real, temp_bias::Real;
        dΔS::Real = LUT_dΔS_DEFAULT)
    F = Sleipnir.Float
    ΔS_max > ΔS_min || throw(ArgumentError(
        "Empty lookup table range: ΔS_max = $(ΔS_max) m is not above ΔS_min = $(ΔS_min) m."))
    dΔS > 0 || throw(ArgumentError("Lookup table resolution dΔS must be positive."))

    n_e = ceil(Int, (ΔS_max - ΔS_min) / dΔS) + 1
    n_w = length(windows)
    PDD = Matrix{F}(undef, n_e, n_w)
    snow = Matrix{F}(undef, n_e, n_w)
    for k in 1:n_w
        window = windows[k]
        @inbounds for i in 1:n_e
            ΔS = F(ΔS_min) + (i - 1) * F(dΔS)
            p, s = pdd_snow_exact(window, ΔS, temp_bias)
            PDD[i, k] = p
            snow[i, k] = s
        end
    end
    # The stored ΔS_max is the last tabulated node, which ceil() may have pushed past the
    # requested one. Recording the actual top of the table keeps the range check exact.
    return ElevationLUT{F}(
        F(ΔS_min), F(ΔS_min) + (n_e - 1) * F(dΔS), F(dΔS), F(1 / dΔS), PDD, snow)
end

@noinline function _lut_range_error(lut::ElevationLUT, ΔS)
    throw(DomainError(ΔS,
        "Elevation offset ΔS = $(ΔS) m is outside the lookup table range " *
        "[$(lut.ΔS_min), $(lut.ΔS_max)] m (glacier elevations ± $(LUT_PAD_DEFAULT) m pad). " *
        "Either the surface left that range or the solve diverged; extrapolating would be " *
        "wrong since PDD keeps growing linearly as ΔS decreases."))
end

"""
    lut_lookup(lut::ElevationLUT, ΔS, k::Int)

Positive degree days and solid precipitation at elevation offset `ΔS` in window `k`, by
linear interpolation in the table.

Throws a `DomainError`, rather than clamping, if `ΔS` is outside the tabulated range or not
finite — clamping would quietly return a wrong mass balance.
"""
@inline function lut_lookup(lut::ElevationLUT, ΔS::R, k::Int) where {R <: Real}
    n_e = size(lut.PDD, 1)
    x = (ΔS - lut.ΔS_min) * lut.inv_dΔS
    (isfinite(x) && zero(R) <= x <= R(n_e - 1)) || _lut_range_error(lut, ΔS)
    # The index carries no derivative; sensitivity to ΔS flows through the weight w, which
    # stays of type R.
    i = floor(Int, x)
    i = i > n_e - 2 ? n_e - 2 : i
    w = x - i
    i1 = i + 1
    @inbounds begin
        pdd = lut.PDD[i1, k] + w * (lut.PDD[i1 + 1, k] - lut.PDD[i1, k])
        snow = lut.snow[i1, k] + w * (lut.snow[i1 + 1, k] - lut.snow[i1, k])
    end
    return pdd, snow
end

"""
    lut_lookup_∂(lut::ElevationLUT, ΔS, k)

Same as [`lut_lookup`](@ref) but also returns the derivatives of `PDD` and `snow` with
respect to `ΔS`. The table is piecewise linear, so the slope of the bracketing segment is
the exact derivative, not an approximation of it.
"""
@inline function lut_lookup_∂(lut::ElevationLUT, ΔS::R, k::Int) where {R <: Real}
    n_e = size(lut.PDD, 1)
    x = (ΔS - lut.ΔS_min) * lut.inv_dΔS
    (isfinite(x) && zero(R) <= x <= R(n_e - 1)) || _lut_range_error(lut, ΔS)
    i = floor(Int, x)
    i = i > n_e - 2 ? n_e - 2 : i
    w = x - i
    i1 = i + 1
    @inbounds begin
        ΔPDD = lut.PDD[i1 + 1, k] - lut.PDD[i1, k]
        Δsnow = lut.snow[i1 + 1, k] - lut.snow[i1, k]
        pdd = lut.PDD[i1, k] + w * ΔPDD
        snow = lut.snow[i1, k] + w * Δsnow
    end
    return pdd, snow, ΔPDD * lut.inv_dΔS, Δsnow * lut.inv_dΔS
end

"""
    MBcache{F <: AbstractFloat}

Everything a mass balance model needs to be evaluated as a source term inside the ice flow
RHS, precomputed once per glacier at `init_cache` time.

The RHS is called many times per mass balance window and must not touch `Rasters`, so the
climate is sliced up front into [`ClimateWindow`](@ref)s and, for models whose rate depends
on the surface only through `ΔS`, collapsed further into an [`ElevationLUT`](@ref).

An empty cache (no windows, empty table) is built when mass balance is off, keeping
`ModelCache` concretely typed either way — a run without mass balance pays neither the
memory nor the build time. Use [`mb_cache_active`](@ref) to tell them apart.

# Fields

  - `windows::Vector{ClimateWindow{F}}`: Climate per mass balance window, in time order.
  - `lut::ElevationLUT{F}`: Elevation lookup table, empty for `:general` models.
  - `ṁ::Matrix{F}`: Buffer holding the most recently computed rate, in m/yr.
  - `t₀::F`, `step_MB::F`: Start of the run and the mass balance step, in decimal years.
  - `ref_hgt::F`: Reference elevation of the climate data, in m.
  - `temp_bias::F`: Temperature bias baked into `lut`, in °C.
  - `H_acc::F`, `H_acc_w::F`: Ice thickness at which accumulation switches on, and the width
    of that switch, in m. The ramp is centred on `H_acc`: it is zero below `H_acc - H_acc_w/2`
    and full above `H_acc + H_acc_w/2`.
  - `H_abl::F`: Ramp width for ablation, in m, rising from zero at a bare bed.
"""
struct MBcache{F <: AbstractFloat}
    windows::Vector{ClimateWindow{F}}
    lut::ElevationLUT{F}
    ṁ::Matrix{F}
    t₀::F
    step_MB::F
    ref_hgt::F
    temp_bias::F
    H_acc::F
    H_acc_w::F
    H_abl::F
    ∂ṁ_max::F
end

"""
    mb_cache_active(cache)

Whether `cache` actually carries precomputed climate, i.e. whether mass balance is being
evaluated inside the ice flow RHS for this simulation. `nothing`, the cache of a mass
balance model that has no RHS form, is never active.
"""
mb_cache_active(cache::MBcache) = !isempty(cache.windows)
mb_cache_active(::Nothing) = false

"""
    lut_∂ṁ_bound(lut::ElevationLUT, mb_model, step_MB, H_abl, H_acc_w)

Upper bound on `|∂ṁ/∂H|` over the whole table, in yr⁻¹.

This is the mass balance contribution to the spectral radius of the ice flow right hand side.
It's a state-independent bound, not the exact maximum at the current state: a stabilised
solver only needs an upper bound to size its stages, and the adjoint hands the RHS an
augmented `[H; θ]` vector that a state-dependent version couldn't read anyway.

From `ṁ = rate(ΔS)·ramp(H)`, with `ramp ≤ 1` and `ramp' ≤ 1.5/w` for a ramp of width `w`,
taking the narrower of the accumulation and ablation ramps:

```math
|∂ṁ/∂H| ≤ \\max|∂rate/∂ΔS| + \\max|rate| · 1.5 / \\min(H_{abl}, H_{acc,w})
```

Both maxima come straight from the table, whose rows are the tabulated elevations and whose
slopes between them are exact.
"""
function lut_∂ṁ_bound(
        lut::ElevationLUT{F}, mb_model, step_MB::F, H_abl::F, H_acc_w::F) where {F}
    isempty(lut) && return zero(F)
    DDF = F(mb_model.DDF)
    c_prcp = F(PRECIP_UNIT_CONVERSION * mb_model.prcp_fac)
    inv_step = one(F) / step_MB

    max_rate = zero(F)
    max_∂rate = zero(F)
    n_e, n_w = size(lut.PDD)
    @inbounds for k in 1:n_w, i in 1:n_e

        r = abs(c_prcp * lut.snow[i, k] - DDF * lut.PDD[i, k]) * inv_step
        max_rate = r > max_rate ? r : max_rate
        if i < n_e
            ∂r = abs(c_prcp * (lut.snow[i + 1, k] - lut.snow[i, k]) -
                     DDF * (lut.PDD[i + 1, k] - lut.PDD[i, k])) * lut.inv_dΔS * inv_step
            max_∂rate = ∂r > max_∂rate ? ∂r : max_∂rate
        end
    end
    # Whichever ramp is narrower sets the steepest `∂ramp/∂H`
    return max_∂rate + max_rate * F(1.5) / min(H_abl, H_acc_w)
end

function _empty_mb_cache(F::Type{<:AbstractFloat} = Sleipnir.Float)
    lut = ElevationLUT{F}(zero(F), zero(F), one(F), one(F),
        Matrix{F}(undef, 0, 0), Matrix{F}(undef, 0, 0))
    return MBcache{F}(ClimateWindow{F}[], lut, Matrix{F}(undef, 0, 0),
        zero(F), one(F), zero(F), zero(F),
        H_ACC_DEFAULT, H_ACC_W_DEFAULT, H_ABL_DEFAULT, zero(F))
end

"""
    window_index(cache::MBcache, t)

Index of the mass balance window containing time `t`.

Window `k` spans `(t₀ + (k-1)·step, t₀ + k·step]`. Times at or before `t₀` map to the first
window and times past the end map to the last.
"""
@inline function window_index(cache::MBcache, t::Real)
    # The nudge makes the interval right-closed in floating point. Landing on a window edge
    # is common, not rare — those times are solver tstops, since the source term is
    # discontinuous in t there. Without it, which window t_k belongs to comes down to rounding.
    k = ceil(Int, (t - cache.t₀) / cache.step_MB - 1e-9)
    n = length(cache.windows)
    return k < 1 ? 1 : (k > n ? n : k)
end

"""
    MB_rate!(ṁ, H, cache::MBcache, mb_model::TImodel1, glacier, t)

Write the instantaneous mass balance rate `ṁ(H, t)`, in m/yr, into `ṁ`.

Rate form of `compute_MB`: that returns the balance accumulated over a window, for
calibration and diagnostics; this returns the rate at which it accrues, to add to `∂H/∂t`.
Integrating the rate over a window with `H` frozen reproduces `compute_MB` exactly; letting
`H` evolve within the window is the operator-splitting error this avoids.

`ṁ` depends on `H` through `S = B + H`, so `∂ṁ/∂H` (the elevation feedback) reaches the
gradient with respect to ice flow parameters like `A` even when no mass balance parameter is
trained.

`DDF` and `prcp_fac` come from `mb_model` on every call, so they can change without
invalidating the lookup table; `temp_bias` is baked into the table and checked against the
cache instead.

`ṁ` is overwritten, not accumulated into — the caller combines it with the rest of the RHS.
"""
function MB_rate!(
        ṁ, H, cache::MBcache, mb_model::TImodel1,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    mb_cache_active(cache) || throw(ArgumentError(
        "MB_rate! called with an inactive mass balance cache. The cache is only " *
        "populated when use_MB is true."))
    get_temp_bias(mb_model) == cache.temp_bias || throw(ArgumentError(
        "Mass balance model temp_bias = $(get_temp_bias(mb_model)) °C does not match " *
        "the $(cache.temp_bias) °C baked into the lookup table. Rebuild the cache."))

    k = window_index(cache, t)
    B = glacier.B
    DDF = mb_model.DDF
    # compute_MB divides the window total by step/(1/12) to express it per month; the rate
    # form instead divides by the window length to express it per year.
    inv_step = 1 / cache.step_MB
    c_prcp = PRECIP_UNIT_CONVERSION * mb_model.prcp_fac
    lut = cache.lut
    ref_hgt = cache.ref_hgt
    H_acc = cache.H_acc
    inv_H_acc_w = 1 / cache.H_acc_w
    inv_H_abl = 1 / cache.H_abl

    @inbounds for j in axes(H, 2), i in axes(H, 1)

        h = H[i, j]
        h = h > 0 ? h : zero(h)
        ΔS = B[i, j] + h - ref_hgt
        pdd, snow = lut_lookup(lut, ΔS, k)
        rate = (c_prcp * snow - DDF * pdd) * inv_step
        # Accumulation switches on around H_acc, ablation rises from a bare bed
        x = rate > 0 ? (h - H_acc) * inv_H_acc_w + oftype(h, 0.5) : h * inv_H_abl
        ṁ[i, j] = rate * smoothstep(x)
    end
    return nothing
end

"""
    MB_rate_∂H!(∂ṁ, H, cache::MBcache, mb_model, glacier, t)

Fill `∂ṁ` with the derivative of the mass balance rate with respect to the ice thickness,
in-place.

`ṁ` at a cell depends only on that cell's `H`, through `S = B + H` and through the ramp, so
the Jacobian is diagonal and a vector-Jacobian product is just an elementwise multiply.

From `ṁ = rate(ΔS) · ramp(H)` and `∂ΔS/∂H = 1`, the product rule gives

```math
∂ṁ/∂H = \\frac{∂rate}{∂ΔS} ramp + rate \\frac{∂ramp}{∂H}
```

Both terms vanish below the ice margin, so the derivative is continuous there. Only manual
adjoints need this; the automatic sensitivity path differentiates [`MB_rate!`](@ref)
directly.
"""
function MB_rate_∂H!(
        ∂ṁ, H, cache::MBcache, mb_model::TImodel1,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    mb_cache_active(cache) || throw(ArgumentError(
        "MB_rate_∂H! called with an inactive mass balance cache. The cache is only " *
        "populated when use_MB is true."))
    get_temp_bias(mb_model) == cache.temp_bias || throw(ArgumentError(
        "Mass balance model temp_bias = $(get_temp_bias(mb_model)) °C does not match " *
        "the $(cache.temp_bias) °C baked into the lookup table. Rebuild the cache."))

    k = window_index(cache, t)
    B = glacier.B
    DDF = mb_model.DDF
    inv_step = 1 / cache.step_MB
    c_prcp = PRECIP_UNIT_CONVERSION * mb_model.prcp_fac
    lut = cache.lut
    ref_hgt = cache.ref_hgt
    H_acc = cache.H_acc
    inv_H_acc_w = 1 / cache.H_acc_w
    inv_H_abl = 1 / cache.H_abl

    @inbounds for j in axes(H, 2), i in axes(H, 1)

        h = H[i, j]
        if h <= 0
            ∂ṁ[i, j] = zero(eltype(∂ṁ))
            continue
        end
        ΔS = B[i, j] + h - ref_hgt
        pdd, snow, ∂pdd, ∂snow = lut_lookup_∂(lut, ΔS, k)
        rate = (c_prcp * snow - DDF * pdd) * inv_step
        ∂rate = (c_prcp * ∂snow - DDF * ∂pdd) * inv_step
        inv_H = rate > 0 ? inv_H_acc_w : inv_H_abl
        x = rate > 0 ? (h - H_acc) * inv_H_acc_w + oftype(h, 0.5) : h * inv_H_abl
        ∂ṁ[i, j] = ∂rate * smoothstep(x) + rate * smoothstep_∂(x) * inv_H
    end
    return nothing
end
function MB_rate_∂H!(∂ṁ, H, cache::MBcache, mb_model::MBmodel,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    throw(ArgumentError(
        "Mass balance model $(typeof(mb_model)) has no ice thickness derivative for the " *
        "ice flow RHS (mb_S_dependence = :$(mb_S_dependence(mb_model)))."))
end

"""
    MB_rate_∂H_maxabs(H, cache::MBcache, mb_model, glacier, t)

Largest `|∂ṁ/∂H|` over the grid, in yr⁻¹, without materialising the derivative.

Feeds the mass balance contribution to the spectral radius, for a stabilised solver to size
its stages. A reduction rather than [`MB_rate_∂H!`](@ref) followed by `maximum`, because it
runs inside a differentiated region every step — allocating and writing a buffer there is
what upsets reverse mode.
"""
function MB_rate_∂H_maxabs(H, cache::MBcache, mb_model::TImodel1,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    mb_cache_active(cache) || return zero(eltype(H))

    k = window_index(cache, t)
    B = glacier.B
    DDF = mb_model.DDF
    inv_step = 1 / cache.step_MB
    c_prcp = PRECIP_UNIT_CONVERSION * mb_model.prcp_fac
    lut = cache.lut
    ref_hgt = cache.ref_hgt
    H_acc = cache.H_acc
    inv_H_acc_w = 1 / cache.H_acc_w
    inv_H_abl = 1 / cache.H_abl

    length(H) == length(B) || throw(DimensionMismatch(
        "Ice thickness has $(length(H)) entries but the bed has $(length(B))."))

    λ = zero(eltype(H))
    # Linear indexing, not (i, j): the solver hands back a flat vector, and the bed shares
    # its layout. 2D indexing here reads out of bounds, which `@inbounds` turns into silent
    # garbage rather than an error.
    @inbounds for idx in eachindex(H)
        h = H[idx]
        h > 0 || continue
        ΔS = B[idx] + h - ref_hgt
        pdd, snow, ∂pdd, ∂snow = lut_lookup_∂(lut, ΔS, k)
        rate = (c_prcp * snow - DDF * pdd) * inv_step
        ∂rate = (c_prcp * ∂snow - DDF * ∂pdd) * inv_step
        inv_H = rate > 0 ? inv_H_acc_w : inv_H_abl
        x = rate > 0 ? (h - H_acc) * inv_H_acc_w + oftype(h, 0.5) : h * inv_H_abl
        a = abs(∂rate * smoothstep(x) + rate * smoothstep_∂(x) * inv_H)
        λ = a > λ ? a : λ
    end
    return λ
end
MB_rate_∂H_maxabs(H, cache::MBcache, mb_model::MBmodel, glacier, t::Real) = zero(eltype(H))

"""
    MB_rate!(ṁ, H, cache::MBcache, mb_model::MBmodel, glacier, t)

Fallback for mass balance models that have no RHS form yet.
"""
function MB_rate!(ṁ, H, cache::MBcache, mb_model::MBmodel,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    throw(ArgumentError(
        "Mass balance model $(typeof(mb_model)) cannot be evaluated as a source term in " *
        "the ice flow RHS (mb_S_dependence = :$(mb_S_dependence(mb_model))). Implement " *
        "MB_rate! for this model."))
end

"""
    init_mb_cache(mb_model::TImodel1, simulation, glacier_idx, θ)

Build the [`MBcache`](@ref) for one glacier. Returns an empty cache when mass balance is off.

The lookup table spans the elevations the glacier can reach — lowest bed to highest initial
surface — padded either side, to absorb thickening beyond the initial state or a bed
emerging below the current minimum. Going outside that range throws rather than
extrapolating.
"""
function init_mb_cache(mb_model::TImodel1, simulation, glacier_idx::Integer, θ)
    F = Sleipnir.Float
    # A stand-in simulation with no parameters isn't requesting the continuous scheme
    hasproperty(simulation, :parameters) || return _empty_mb_cache(F)
    simparams = simulation.parameters.simulation
    simparams.use_MB || return _empty_mb_cache(F)

    glacier = simulation.glaciers[glacier_idx]
    tspan = simparams.tspan
    step_MB = simparams.step_MB
    temp_bias = F(get_temp_bias(mb_model))

    windows = precompute_climate_windows(glacier.climate, tspan, step_MB)
    ref_hgt = F(windows[begin].step.ref_hgt)

    # Reachable elevation offsets: from a fully deglaciated bed to the highest surface the
    # glacier starts with, plus a pad either side.
    ΔS_min = F(minimum(glacier.B)) - ref_hgt - LUT_PAD_DEFAULT
    ΔS_max = F(maximum(glacier.B .+ glacier.H₀)) - ref_hgt + LUT_PAD_DEFAULT
    lut = build_elevation_lut(windows, ΔS_min, ΔS_max, temp_bias)

    return MBcache{F}(
        windows, lut, zeros(F, glacier.nx, glacier.ny),
        F(tspan[1]), F(step_MB), ref_hgt, temp_bias,
        H_ACC_DEFAULT, H_ACC_W_DEFAULT, H_ABL_DEFAULT,
        lut_∂ṁ_bound(lut, mb_model, F(step_MB), H_ABL_DEFAULT, H_ACC_W_DEFAULT))
end

mb_cache_type(::TImodel1) = MBcache{Sleipnir.Float}

# Display setup
Base.show(io::IO, ::MIME"text/plain", cache::MBcache) = Base.show(io, cache)
function Base.show(io::IO, cache::MBcache)
    if !mb_cache_active(cache)
        print(io, "MBcache (inactive: mass balance is not evaluated in the ice flow RHS)")
        return
    end
    println(io, "MBcache")
    println(io, "   windows   = ", length(cache.windows),
        " × ", cache.step_MB, " yr from ", cache.t₀)
    println(io, "   LUT       = ", size(cache.lut.PDD, 1), " elevations over [",
        cache.lut.ΔS_min, ", ", cache.lut.ΔS_max, "] m at ", cache.lut.dΔS, " m")
    println(io, "   ref_hgt   = ", cache.ref_hgt, " m")
    println(io, "   temp_bias = ", cache.temp_bias, " °C")
    print(io, "   ramps     = H_acc ", cache.H_acc, " ± ", cache.H_acc_w / 2,
        " m, H_abl ", cache.H_abl, " m")
end
