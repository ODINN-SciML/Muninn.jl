
"""
    plot_cumulative_mb(
        mb_model::MBmodel,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0),
        kwargs...,
    ) where {F <: AbstractFloat}

Plot the mean annual mass balance map (m w.e. yr⁻¹) for a single glacier over
`[t_start, t_end]` using a calibrated mass balance model and static geometry.

Extends [`Sleipnir.plot_cumulative_mb`](@ref) with a dispatch that accepts an
[`MBmodel`](@ref) and a glacier directly.  The cumulative MB is computed via
[`compute_cumulative_MB`](@ref) and divided by `(t_end - t_start)` to obtain an
annual rate before plotting.  All `kwargs` are forwarded (e.g. `title`, `colormap`).
"""
function Sleipnir.plot_cumulative_mb(
        mb_model::MBmodel,
        glacier::AbstractGlacier,
        t_start::F,
        t_end::F;
        step::F = F(1.0 / 12.0),
        kwargs...) where {F <: AbstractFloat}
    cumulative_mb = compute_cumulative_MB(mb_model, glacier, t_start, t_end; step)
    n_years = t_end - t_start
    annual_mb = cumulative_mb ./ n_years
    # Build a minimal Results from the glacier for the Sleipnir plotting machinery.
    # MBmodel <: AbstractModel satisfies the constructor's ifm type constraint;
    # S is provided explicitly so ifm.S is never accessed.
    results = Sleipnir.Results(glacier, mb_model;
        S = glacier.S,
        H = [glacier.H₀],
        MB = [annual_mb],
        tspan = (Sleipnir.Float(t_start), Sleipnir.Float(t_end)))
    return Sleipnir.plot_cumulative_mb(results; colorbar_label = "m w.e. yr⁻¹", kwargs...)
end
