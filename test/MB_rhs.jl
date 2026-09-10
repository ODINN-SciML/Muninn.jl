# Tests for mass balance as a source term in the ice flow RHS: the cache, the elevation
# lookup table, and the rate itself.

# init_mb_cache only reads `.glaciers` and `.parameters` off the simulation, so a stand-in
# keeps these tests in Muninn instead of pulling in a Huginn Prediction.
struct fakeSimulation{G, P}
    glaciers::Vector{G}
    parameters::P
end

function _mb_rhs_setup(; MB_scheme::Symbol = :continuous, tspan = (2010.0, 2011.0),
        step_MB = 1.0/12.0, temp_bias = 0.0)
    rgi_ids = ["RGI60-11.03638"]
    rgi_paths = get_rgi_paths()
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(simulation = SimulationParameters(
        use_MB = true, MB_scheme = MB_scheme, use_velocities = false,
        tspan = tspan, step_MB = step_MB, test_mode = true,
        multiprocessing = false, rgi_paths = rgi_paths))
    glacier = initialize_glaciers(rgi_ids, params)[1]
    glacier.S .= glacier.B .+ glacier.H₀
    mb_model = TImodel1(params; temp_bias = temp_bias)
    simulation = fakeSimulation([glacier], params)
    return (; params, glacier, mb_model, simulation, tspan, step_MB)
end

function smoothstep_test()
    # Endpoints and clamping.
    @test smoothstep(0.0) == 0.0
    @test smoothstep(1.0) == 1.0
    @test smoothstep(-3.0) == 0.0
    @test smoothstep(7.0) == 1.0
    @test smoothstep(0.5) == 0.5

    # Monotone on the ramp.
    xs = range(0.0, 1.0; length = 101)
    @test issorted(smoothstep.(xs))

    # C¹ at both ends is the whole point: a hard mask is a step discontinuity that an
    # adaptive controller cannot resolve and an adjoint cannot differentiate.
    h = 1e-6
    @test abs((smoothstep(h) - smoothstep(0.0)) / h) < 1e-5
    @test abs((smoothstep(1.0) - smoothstep(1.0 - h)) / h) < 1e-5

    # Derivative is bounded by 3/2 over the ramp, as documented.
    d = [(smoothstep(x + h) - smoothstep(x)) / h for x in range(0.0, 1.0 - h; length = 501)]
    @test maximum(d) <= 1.5 + 1e-4

    # Type generic, so the same code path serves AD.
    @test smoothstep(0.5f0) isa Float32
end

function mb_S_dependence_test()
    params = Parameters(simulation = SimulationParameters(multiprocessing = false))
    @test mb_S_dependence(TImodel1(params)) == :elevation_only
    # Anything that has not opted in must fall back to the conservative default.
    @test mb_S_dependence(TImodel2(params)) == :general
    @test mb_S_dependence(DummyMBDefault()) == :general
end

"""
Check the lookup table against the exact daily sums it approximates.

`PDD` and `snow` are exactly piecewise linear in `ΔS`, so linear interpolation is exact
except within one table cell of a breakpoint. This bounds the residual error rather than
assuming it.
"""
function elevation_lut_test()
    s = _mb_rhs_setup()
    windows = precompute_climate_windows(s.glacier.climate, s.tspan, s.step_MB)
    lut = build_elevation_lut(windows, -1500.0, 1500.0, 0.0)

    @test size(lut.PDD) == size(lut.snow)
    @test size(lut.PDD, 2) == length(windows)
    @test lut.ΔS_min == -1500.0
    @test lut.ΔS_max >= 1500.0
    @test !isempty(lut)

    # Table nodes are exact by construction.
    for k in (1, length(windows))
        for ΔS in (-1500.0, -700.0, 0.0, 1200.0)
            pdd_ref, snow_ref = Muninn.pdd_snow_exact(windows[k], ΔS, 0.0)
            pdd, snow = Muninn.lut_lookup(lut, ΔS, k)
            @test pdd ≈ pdd_ref rtol=1e-12
            @test snow ≈ snow_ref rtol=1e-12
        end
    end

    # Between nodes, the interpolation error must stay negligible on the quantity that
    # actually matters: the monthly mass balance in m.
    mb_model = s.mb_model
    max_err = 0.0
    for k in 1:length(windows), ΔS in range(-1490.0, 1490.0; length = 401)

        pdd_ref, snow_ref = Muninn.pdd_snow_exact(windows[k], ΔS, 0.0)
        pdd, snow = Muninn.lut_lookup(lut, ΔS, k)
        mb_ref = PRECIP_UNIT_CONVERSION * mb_model.prcp_fac * snow_ref -
                 mb_model.DDF * pdd_ref
        mb = PRECIP_UNIT_CONVERSION * mb_model.prcp_fac * snow - mb_model.DDF * pdd
        max_err = max(max_err, abs(mb - mb_ref))
    end
    @test max_err < 1e-3

    # Out of range throws instead of extrapolating: constant extrapolation is exact only
    # above the range, since PDD keeps growing linearly as ΔS decreases.
    @test_throws DomainError Muninn.lut_lookup(lut, -1500.1, 1)
    @test_throws DomainError Muninn.lut_lookup(lut, lut.ΔS_max + 0.1, 1)
    # A diverging solve sends ΔS to absurd values; it must surface here, not as an
    # InexactError inside the index conversion.
    @test_throws DomainError Muninn.lut_lookup(lut, 1e24, 1)
    @test_throws DomainError Muninn.lut_lookup(lut, NaN, 1)

    @test_throws ArgumentError build_elevation_lut(windows, 100.0, 100.0, 0.0)
    @test_throws ArgumentError build_elevation_lut(windows, -10.0, 10.0, 0.0; dΔS = 0.0)
end

"""
Check that the cache declares the type it returns, and that models without an RHS form are
refused rather than silently mis-evaluated.

`cache_type(model::Model)` becomes a type parameter of `Prediction` and `Inversion`, so a
mismatch between `mb_cache_type` and `init_mb_cache` is a type instability in the RHS hot
loop rather than an error anywhere obvious.
"""
function mb_cache_init_test()
    s = _mb_rhs_setup(MB_scheme = :continuous)
    cache = Sleipnir.init_mb_cache(s.mb_model, s.simulation, 1, nothing)
    @test cache isa MBcache
    @test mb_cache_active(cache)
    @test typeof(cache) == Sleipnir.mb_cache_type(s.mb_model)
    @test length(cache.windows) == 12
    @test cache.t₀ == s.tspan[1]
    @test cache.step_MB == s.step_MB
    @test size(cache.ṁ) == (s.glacier.nx, s.glacier.ny)
    @test cache.temp_bias == 0.0

    # The table must bracket every elevation offset the glacier can reach, with the pad.
    ΔS_lo = minimum(s.glacier.B) - cache.ref_hgt
    ΔS_hi = maximum(s.glacier.B .+ s.glacier.H₀) - cache.ref_hgt
    @test cache.lut.ΔS_min <= ΔS_lo - Muninn.LUT_PAD_DEFAULT + 1e-6
    @test cache.lut.ΔS_max >= ΔS_hi + Muninn.LUT_PAD_DEFAULT - 1e-6

    # Models with no RHS form are refused explicitly.
    @test_throws ArgumentError MB_rate!(
        cache.ṁ, s.glacier.H₀, cache, TImodel2(s.params), s.glacier, 2010.5)

    # A model whose temp_bias no longer matches the one baked into the table is refused.
    stale = TImodel1(s.params; temp_bias = 1.5)
    @test_throws ArgumentError MB_rate!(
        cache.ṁ, s.glacier.H₀, cache, stale, s.glacier, 2010.5)

    # Scoped to the MB_scheme staging flag: while both schemes coexist, the discrete one
    # must pay neither the memory nor the build time, yet keep the cache type identical so
    # ModelCache stays concretely typed. Delete with the flag.
    s_disc = _mb_rhs_setup(MB_scheme = :discrete)
    cache_disc = Sleipnir.init_mb_cache(s_disc.mb_model, s_disc.simulation, 1, nothing)
    @test typeof(cache_disc) == typeof(cache)
    @test !mb_cache_active(cache_disc)
    @test isempty(cache_disc.lut)
    @test_throws ArgumentError MB_rate!(
        cache_disc.ṁ, s_disc.glacier.H₀, cache_disc, s_disc.mb_model,
        s_disc.glacier, 2010.5)
end

function mb_window_index_test()
    s = _mb_rhs_setup()
    cache = Sleipnir.init_mb_cache(s.mb_model, s.simulation, 1, nothing)
    t₀, step = s.tspan[1], s.step_MB
    n = length(cache.windows)

    # Window k is (t₀ + (k-1)·step, t₀ + k·step]: the right edge belongs to window k.
    for k in 1:n
        @test Muninn.window_index(cache, t₀ + k * step) == k
        @test Muninn.window_index(cache, t₀ + (k - 0.5) * step) == k
    end
    # Edges are solver tstops, so landing exactly on one is the common case and must not
    # depend on the last bit of the division.
    for k in 1:(n - 1)
        @test Muninn.window_index(cache, t₀ + k * step + 1e-12) == k
    end
    # Outside the span, clamp to the ends rather than throwing.
    @test Muninn.window_index(cache, t₀) == 1
    @test Muninn.window_index(cache, t₀ - 1.0) == 1
    @test Muninn.window_index(cache, t₀ + (n + 5) * step) == n
end

"""
The load-bearing test: integrating the mass balance *rate* over a window with a frozen
surface must reproduce the window total that `compute_MB` returns.

This is not a discrete-versus-continuous transition check, and it does not expire with the
`MB_scheme` flag. `compute_MB` and `downscale_2D_climate` are the kernel Muninn calibrates
`DDF`, `prcp_fac` and `temp_bias` against geodetic observations with (`calibration.jl`),
integrating mass balance over the Hugonnet window with a fixed surface and no ice flow.
If `MB_rate!` ever drifts from them, the model being calibrated stops being the model being
integrated — a silent scientific error with no other guard. Freezing the surface here
matches how the calibration uses them, and isolates the formulation from the
operator-splitting difference, which only appears once `H` evolves within a window.

The comparison is restricted to cells where both ramps are saturated. Where they are not,
the two are deliberately different: `apply_MB_mask!` applies a hard mask, `MB_rate!` a C¹
ramp.
"""
function mb_rate_matches_compute_MB_test()
    s = _mb_rhs_setup()
    glacier, mb_model, step_MB = s.glacier, s.mb_model, s.step_MB
    cache = Sleipnir.init_mb_cache(mb_model, s.simulation, 1, nothing)

    H = glacier.H₀
    glacier.S .= glacier.B .+ H
    # Saturated ramps, so ramp == 1 for accumulation and ablation alike. The accumulation
    # ramp is centred on H_acc, so it only saturates half a width above it.
    saturated = H .>= Muninn.H_ACC_DEFAULT + Muninn.H_ACC_W_DEFAULT / 2
    @test count(saturated) > 100

    for k in (1, 6, 12)
        t_k = s.tspan[1] + k * step_MB

        # The window total, computed exactly as calibration.jl computes it.
        get_cumulative_climate!(glacier.climate, Sleipnir.Float(t_k),
            Sleipnir.Float(step_MB))
        downscale_2D_climate!(glacier; temp_bias = get_temp_bias(mb_model))
        MB_total = compute_MB(mb_model, glacier.climate.climate_2D_step, step_MB)

        # The same quantity as a rate, integrated over the window.
        MB_rate!(cache.ṁ, H, cache, mb_model, glacier, t_k)
        MB_integrated = cache.ṁ .* step_MB

        Δ = abs.(MB_integrated[saturated] .- MB_total[saturated])
        @test maximum(Δ) < 1e-3
    end
end

"""
Check the ramps, and the elevation feedback that survives once they saturate.

The feedback is the reason this whole formulation needs an adjoint at all: no mass balance
parameter is trained, but `ṁ` depends on `H` through `S = B + H`, so `∂ṁ/∂H` sits in the
adjoint Jacobian and reaches the gradient with respect to ice flow parameters such as `A`.
"""
function mb_rate_ramp_test()
    s = _mb_rhs_setup()
    glacier, mb_model = s.glacier, s.mb_model
    cache = Sleipnir.init_mb_cache(mb_model, s.simulation, 1, nothing)
    t = s.tspan[1] + 6 * s.step_MB
    sz = size(glacier.H₀)

    rate_at(h, i, j) = begin
        MB_rate!(cache.ṁ, fill(h, sz), cache, mb_model, glacier, t)
        cache.ṁ[i, j]
    end

    # A bare bed produces no mass balance at all, whatever the climate says.
    MB_rate!(cache.ṁ, zeros(sz), cache, mb_model, glacier, t)
    @test all(iszero, cache.ṁ)

    # Negative H is treated as bare ice, not extrapolated below the bed.
    ṁ_zero = copy(cache.ṁ)
    MB_rate!(cache.ṁ, fill(-5.0, sz), cache, mb_model, glacier, t)
    @test cache.ṁ == ṁ_zero

    # Pick an ablating and an accumulating cell at a thickness where both ramps saturate,
    # so the two branches of the ramp are each exercised on a cell that uses them.
    MB_rate!(cache.ṁ, fill(Muninn.H_ACC_DEFAULT, sz), cache, mb_model, glacier, t)
    i_abl, j_abl = Tuple(argmin(cache.ṁ))
    i_acc, j_acc = Tuple(argmax(cache.ṁ))
    @test cache.ṁ[i_abl, j_abl] < 0

    # Over the ablation ramp the magnitude grows smoothly from zero to its full value.
    ramp_hs = range(0.0, Muninn.H_ABL_DEFAULT; length = 9)
    ramped = [rate_at(h, i_abl, j_abl) for h in ramp_hs]
    @test ramped[1] == 0.0
    @test issorted(abs.(ramped))
    @test all(isfinite, ramped)

    # Past saturation the ramp is constant, so what remains is the elevation feedback: a
    # thicker glacier has a higher, colder surface, so ṁ is non-decreasing in H. This holds
    # for both branches — less melt where it ablates, more snow where it accumulates.
    saturated_hs = (Muninn.H_ACC_DEFAULT, 25.0, 50.0, 100.0, 200.0)
    for (i, j) in ((i_abl, j_abl), (i_acc, j_acc))
        feedback = [rate_at(h, i, j) for h in saturated_hs]
        @test all(isfinite, feedback)
        @test issorted(feedback)
        # And it is a real dependence, not a rounding artefact.
        @test feedback[end] - feedback[1] > 1e-6
    end

    # The buffer is overwritten, never accumulated into, so the caller owns how the source
    # term is combined with the rest of the RHS.
    MB_rate!(cache.ṁ, glacier.H₀, cache, mb_model, glacier, t)
    first_pass = copy(cache.ṁ)
    MB_rate!(cache.ṁ, glacier.H₀, cache, mb_model, glacier, t)
    @test cache.ṁ == first_pass
end

function mb_rate_∂H_test()
    s = _mb_rhs_setup()
    glacier, mb_model = s.glacier, s.mb_model
    cache = Sleipnir.init_mb_cache(mb_model, s.simulation, 1, nothing)
    t = s.tspan[1] + 6 * s.step_MB
    sz = size(glacier.H₀)
    lut = cache.lut

    ∂ = similar(glacier.H₀)
    ṁp = similar(glacier.H₀)
    ṁm = similar(glacier.H₀)

    # Cells sitting on a lookup table knot are excluded throughout: the table is piecewise
    # linear there, so a two-sided difference straddles two segments and measures neither
    # one-sided slope. That is a genuine kink, not an error in the derivative.
    lut_interior(H) = begin
        x = (glacier.B .+ H .- cache.ref_hgt .- lut.ΔS_min) .* lut.inv_dΔS
        0.1 .< (x .- floor.(x)) .< 0.9
    end

    # Past both ramps smoothstep saturates, so ṁ is exactly linear in H inside a segment
    # and the central difference is exact to roundoff.
    H = fill(120.0, sz)
    MB_rate_∂H!(∂, H, cache, mb_model, glacier, t)
    keep = lut_interior(H)
    # An absolute count, not a fraction: B is clustered modulo the table spacing, so a
    # uniform shift in H moves every cell's position within its segment together.
    @test count(keep) > 1000

    ε = 1e-3
    MB_rate!(ṁp, H .+ ε, cache, mb_model, glacier, t)
    MB_rate!(ṁm, H .- ε, cache, mb_model, glacier, t)
    fd = (ṁp .- ṁm) ./ (2ε)
    # Perturbing the whole field at once also rules out off-diagonal coupling: if ṁ at a
    # cell depended on any neighbour's H, this would not match the diagonal derivative.
    @test maximum(abs.(fd[keep] .- ∂[keep])) < 1e-8
    # And the feedback is a real quantity, not accidentally zero everywhere.
    @test maximum(abs.(∂[keep])) > 1e-4

    # Perturbing one cell isolates one diagonal entry and shows nothing else moves.
    idx = findall(keep)
    for I in idx[round.(Int, range(1, length(idx); length = 5))]
        Hp = copy(H)
        Hp[I] += ε
        Hm = copy(H)
        Hm[I] -= ε
        MB_rate!(ṁp, Hp, cache, mb_model, glacier, t)
        MB_rate!(ṁm, Hm, cache, mb_model, glacier, t)
        @test isapprox((ṁp[I] - ṁm[I]) / (2ε), ∂[I]; atol = 1e-8)
        @test count(!iszero, ṁp .- ṁm) == 1
    end

    # Inside the ramp both factors of rate * ramp depend on H, so this is what actually
    # exercises the product rule.
    H_ramp = fill(0.5 * Muninn.H_ABL_DEFAULT, sz)
    MB_rate_∂H!(∂, H_ramp, cache, mb_model, glacier, t)
    ε_r = 1e-5
    MB_rate!(ṁp, H_ramp .+ ε_r, cache, mb_model, glacier, t)
    MB_rate!(ṁm, H_ramp .- ε_r, cache, mb_model, glacier, t)
    fd_r = (ṁp .- ṁm) ./ (2ε_r)
    # Cells whose rate changes sign across the perturbation swap ramp branch, the other
    # genuine kink, and are excluded for the same reason as the knots.
    keep_r = lut_interior(H_ramp) .& (sign.(ṁp) .== sign.(ṁm))
    @test count(keep_r) > 1000
    @test maximum(abs.(fd_r[keep_r] .- ∂[keep_r])) < 1e-5

    # Below the ice margin the source and its derivative both vanish, so the adjoint sees
    # no feedback on bare ground.
    MB_rate_∂H!(∂, zeros(sz), cache, mb_model, glacier, t)
    @test all(iszero, ∂)
    MB_rate_∂H!(∂, fill(-5.0, sz), cache, mb_model, glacier, t)
    @test all(iszero, ∂)
end
