export MBcache, ElevationLUT, MB_rate!, build_elevation_lut, mb_S_dependence,
       smoothstep, mb_cache_active

import Sleipnir: init_mb_cache, mb_cache_type

# Ice thickness scales (m) over which the mass balance rate is ramped in from zero.
# Accumulation is allowed on a bare bed as soon as a little ice is present; ablation is
# switched off much closer to zero so that a cell cannot melt through the bed.
# Phase 0 measured the final H to be insensitive to H_abl between 0.5 m and 1.0 m
# (RMS 0.009 m) and to move noticeably only at 5.0 m, which fixes the default here.
const H_ACC_DEFAULT = Sleipnir.Float(10.0)
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

Used to fade the mass balance rate in and out with ice thickness. The discrete scheme
applies a hard mask (`apply_MB_mask!`), which is fine for a callback but is a step
discontinuity in `H` when the same quantity becomes a source term in the RHS: an adaptive
error controller cannot resolve it, and an adjoint sees a derivative that is zero almost
everywhere and undefined on a set of measure zero. This ramp is C¹ with a bounded
derivative of at most `1.5` over the ramp width.
"""
@inline function smoothstep(x::R) where {R <: Real}
    x <= zero(R) && return zero(R)
    x >= one(R) && return one(R)
    return x * x * (R(3) - R(2) * x)
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

so a uniform table with linear interpolation is exact except within one cell of a
breakpoint. Phase 0 measured the resulting error on a monthly mass balance at 7.1e-5 m for
`dΔS = 1` m, and the lookup at roughly an eighth of the cost of one `SIA2D!` call, against
ten times that for the equivalent daily loop.

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
        "Elevation offset ΔS = $(ΔS) m falls outside the mass balance lookup table " *
        "range [$(lut.ΔS_min), $(lut.ΔS_max)] m. The table is sized from the elevations " *
        "the glacier can reach plus a $(LUT_PAD_DEFAULT) m pad, so this means either the " *
        "ice surface left that range or the solve diverged. The table is never " *
        "extrapolated silently: constant extrapolation is exact only above the range, " *
        "since PDD keeps growing linearly as ΔS decreases."))
end

"""
    lut_lookup(lut::ElevationLUT, ΔS, k::Int)

Positive degree days and solid precipitation at elevation offset `ΔS` in window `k`, by
linear interpolation in the table.

Throws a `DomainError` if `ΔS` is outside the tabulated range or is not finite, rather than
clamping. Clamping would quietly return a wrong mass balance, and a diverging solve would
then surface much later as an unexplained result.
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
    MBcache{F <: AbstractFloat}

Everything a mass balance model needs to be evaluated as a source term inside the ice flow
RHS, precomputed once per glacier at `init_cache` time.

The RHS is called many times per mass balance window and must not touch `Rasters`, so the
climate is sliced up front into [`ClimateWindow`](@ref)s and, for models whose rate depends
on the surface only through `ΔS`, collapsed further into an [`ElevationLUT`](@ref).

An *empty* cache (no windows, empty table) is what gets built whenever mass balance is not
being evaluated in the RHS — mass balance switched off, or `MB_scheme = :discrete`. Keeping
the type the same either way means `ModelCache` stays concretely typed and the discrete path
pays neither the memory nor the build time. Use [`mb_cache_active`](@ref) to tell them apart.

# Fields

  - `windows::Vector{ClimateWindow{F}}`: Climate per mass balance window, in time order.
  - `lut::ElevationLUT{F}`: Elevation lookup table, empty for `:general` models.
  - `ṁ::Matrix{F}`: Buffer holding the most recently computed rate, in m/yr.
  - `t₀::F`, `step_MB::F`: Start of the run and the mass balance step, in decimal years.
  - `ref_hgt::F`: Reference elevation of the climate data, in m.
  - `temp_bias::F`: Temperature bias baked into `lut`, in °C.
  - `H_acc::F`, `H_abl::F`: Ramp widths for accumulation and ablation, in m.
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
    H_abl::F
end

"""
    mb_cache_active(cache)

Whether `cache` actually carries precomputed climate, i.e. whether mass balance is being
evaluated inside the ice flow RHS for this simulation. `nothing`, the cache of a mass
balance model that has no RHS form, is never active.
"""
mb_cache_active(cache::MBcache) = !isempty(cache.windows)
mb_cache_active(::Nothing) = false

function _empty_mb_cache(F::Type{<:AbstractFloat} = Sleipnir.Float)
    lut = ElevationLUT{F}(zero(F), zero(F), one(F), one(F),
        Matrix{F}(undef, 0, 0), Matrix{F}(undef, 0, 0))
    return MBcache{F}(ClimateWindow{F}[], lut, Matrix{F}(undef, 0, 0),
        zero(F), one(F), zero(F), zero(F), H_ACC_DEFAULT, H_ABL_DEFAULT)
end

"""
    window_index(cache::MBcache, t)

Index of the mass balance window containing time `t`.

Window `k` spans `(t₀ + (k-1)·step, t₀ + k·step]`, matching the times at which the discrete
scheme applies its jumps, so the continuous source term uses exactly the same climate over
exactly the same intervals. Times at or before `t₀` map to the first window and times past
the end map to the last.
"""
@inline function window_index(cache::MBcache, t::Real)
    # The nudge makes the interval right-closed even in floating point. Landing exactly on
    # a window edge is the common case, not a rare one: those times are solver tstops,
    # because the source term is discontinuous in t there. Without it, whether t_k belongs
    # to window k or k+1 would come down to the last bit of (t - t₀)/step.
    k = ceil(Int, (t - cache.t₀) / cache.step_MB - 1e-9)
    n = length(cache.windows)
    return k < 1 ? 1 : (k > n ? n : k)
end

"""
    MB_rate!(ṁ, H, cache::MBcache, mb_model::TImodel1, glacier, t)

Write the instantaneous mass balance rate `ṁ(H, t)`, in m/yr, into `ṁ`.

This is the continuous counterpart of `compute_MB`: where the latter returns the mass
balance accumulated over a whole window, to be applied as a jump, this returns the rate at
which it accrues, to be added to `∂H/∂t`. Integrating the rate over a window with a frozen
surface reproduces `compute_MB` exactly; the two differ once `H` evolves within a window,
which is precisely the operator-splitting error the continuous form removes.

The rate depends on `H` through the surface `S = B + H`, and that dependence is what an
adjoint needs even when no mass balance parameter is being trained: `∂ṁ/∂H` is the elevation
feedback, it sits in the adjoint Jacobian, and so it reaches the gradient with respect to
ice flow parameters such as `A`.

`DDF` and `prcp_fac` are read from `mb_model` on every call rather than from the cache, so
they remain free to change without invalidating the lookup table. `temp_bias` is baked into
the table and is checked against the cache instead.

`ṁ` is fully overwritten, never accumulated into, so the caller owns how it is combined
with the rest of the RHS.
"""
function MB_rate!(
        ṁ, H, cache::MBcache, mb_model::TImodel1,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    mb_cache_active(cache) || throw(ArgumentError(
        "MB_rate! called with an inactive mass balance cache. The cache is only " *
        "populated when use_MB is true and MB_scheme is :continuous."))
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
    inv_H_acc = 1 / cache.H_acc
    inv_H_abl = 1 / cache.H_abl

    @inbounds for j in axes(H, 2), i in axes(H, 1)

        h = H[i, j]
        h = h > 0 ? h : zero(h)
        ΔS = B[i, j] + h - ref_hgt
        pdd, snow = lut_lookup(lut, ΔS, k)
        rate = (c_prcp * snow - DDF * pdd) * inv_step
        ramp = rate > 0 ? smoothstep(h * inv_H_acc) : smoothstep(h * inv_H_abl)
        ṁ[i, j] = rate * ramp
    end
    return nothing
end

"""
    MB_rate!(ṁ, H, cache::MBcache, mb_model::MBmodel, glacier, t)

Fallback for mass balance models that have no RHS form yet.
"""
function MB_rate!(ṁ, H, cache::MBcache, mb_model::MBmodel,
        glacier::Sleipnir.AbstractGlacier, t::Real)
    throw(ArgumentError(
        "Mass balance model $(typeof(mb_model)) cannot be evaluated as a source term in " *
        "the ice flow RHS (mb_S_dependence = :$(mb_S_dependence(mb_model))). Use " *
        "MB_scheme = :discrete for this model."))
end

"""
    init_mb_cache(mb_model::TImodel1, simulation, glacier_idx, θ)

Build the [`MBcache`](@ref) for one glacier.

Returns an empty cache unless mass balance is switched on *and* `MB_scheme` is
`:continuous`, so nothing changes for the discrete path.

The lookup table is sized from the elevations the glacier can reach — the lowest bed and the
highest initial surface — padded either side. The pad is what absorbs a surface that
thickens beyond its initial state or a bed that emerges below the current minimum; going
outside it throws rather than extrapolating.
"""
function init_mb_cache(mb_model::TImodel1, simulation, glacier_idx::Integer, θ)
    F = Sleipnir.Float
    simparams = simulation.parameters.simulation
    (simparams.use_MB && simparams.MB_scheme == :continuous) || return _empty_mb_cache(F)

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
        H_ACC_DEFAULT, H_ABL_DEFAULT)
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
    print(io, "   ramps     = H_acc ", cache.H_acc, " m, H_abl ", cache.H_abl, " m")
end
