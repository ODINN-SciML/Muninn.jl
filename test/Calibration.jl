
function _calibration_test_params(rgi_ids, calibration_tspan)
    rgi_paths = Sleipnir.get_rgi_paths()
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    return Sleipnir.Parameters(
        simulation = Sleipnir.SimulationParameters(
        use_MB = true,
        use_velocities = false,
        tspan = calibration_tspan,
        multiprocessing = false,
        workers = 1,
        test_mode = true,
        rgi_paths = rgi_paths),
    )
end

function calibrate_ti_model_test()
    # Workers spawned by earlier tests only have Sleipnir loaded; they lack
    # Muninn's calibrate_ti_model and would error when pmap dispatches to them.
    # Subsequent tests in this suite use multiprocessing=false so no re-spawn needed.
    nworkers() > 1 && rmprocs(workers())

    rgi_ids = ["RGI60-11.03638"]

    # Calibration period: 2000-2020 (Hugonnet et al.)
    calibration_tspan = (2000.0, 2020.0)

    params = _calibration_test_params(rgi_ids, calibration_tspan)

    glacier = Sleipnir.initialize_glaciers(rgi_ids, params)[1]

    # Attach a synthetic geodetic MB observation to the glacier
    # (–0.5 m w.e. yr⁻¹ is a plausible negative mass balance for this glacier)
    geodetic_mb = Sleipnir.Float(-0.5)
    dhdt_data = Sleipnir.DhdtData(
        (Sleipnir.Float(calibration_tspan[1]), Sleipnir.Float(calibration_tspan[2])),
        geodetic_mb)
    glacier = Sleipnir.Glacier2D(glacier; dhdtData = dhdt_data, geodetic_MB = geodetic_mb)

    # --- compute_mean_annual_MB: verify it returns a finite scalar ---------------
    TI_default = TImodel1(params)
    mb_default = compute_mean_annual_MB(
        TI_default, glacier,
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(mb_default)

    # --- calibrate_ti_model: DDF should drive modelled MB to target --------------
    TI_cal = calibrate_ti_model(glacier, params)
    @test TI_cal isa TImodel1

    # The calibrated DDF must lie within the default search bounds
    @test TI_cal.DDF >= params.physical.DDF_min
    @test TI_cal.DDF <= params.physical.DDF_max

    # prcp_fac should be 1.0 (DDF-only calibration in this test)
    @test TI_cal.prcp_fac ≈ Sleipnir.Float(1.0)

    # temp_bias should be 0.0 (DDF-only calibration in this test)
    @test TI_cal.temp_bias ≈ Sleipnir.Float(0.0)

    # The modelled MB with the calibrated model should be very close to the target
    mb_cal = compute_mean_annual_MB(
        TI_cal, glacier,
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(mb_cal)
    @test abs(mb_cal - geodetic_mb) < 1e-4  # within 0.1 mm w.e. yr⁻¹

    # --- calibrate_MB_model!: high-level entry point expands to a per-glacier vector
    model = Sleipnir.Model(; iceflow = nothing, mass_balance = TImodel1(params))
    calibrate_MB_model!(model, [glacier], params)
    @test model.mass_balance isa Vector{<:TImodel1}
    @test length(model.mass_balance) == 1
    # get_mb_model should return the per-glacier model, matching the direct fit
    @test get_mb_model(model.mass_balance, 1).DDF ≈ TI_cal.DDF

    # --- get_mb_model dispatch: single shared model vs per-glacier vector --------
    single = TImodel1(params)
    @test get_mb_model(single) === single          # single model, any glacier
    @test get_mb_model(single, 7) === single
    @test get_mb_model(model.mass_balance, 1) === model.mass_balance[1]

    # --- calibrate_MB_model! is a no-op without a calibration routine ------------
    # TImodel2 has no calibration method, so the model is returned unchanged.
    ti2 = TImodel2(params)
    model2 = Sleipnir.Model(; iceflow = nothing, mass_balance = ti2)
    calibrate_MB_model!(model2, [glacier], params)
    @test model2.mass_balance === ti2              # left untouched (not vectorized)

    # --- calibration_period override --------------------------------------------
    TI_override = calibrate_ti_model(
        glacier, params;
        calibration_period = (Sleipnir.Float(2005.0), Sleipnir.Float(2015.0)))
    @test TI_override isa TImodel1
    @test isfinite(TI_override.DDF)
end

function calibrate_ti_model_temp_bias_test()
    rgi_ids = ["RGI60-11.03638"]
    calibration_tspan = (2000.0, 2020.0)
    params = _calibration_test_params(rgi_ids, calibration_tspan)
    glacier = Sleipnir.initialize_glaciers(rgi_ids, params)[1]

    # Set a strongly negative target that DDF_max + prcp_fac_min together still
    # can't reach — forcing step 3 (temp_bias). A large positive temp_bias (warmer
    # climate) increases melt and converts snow to rain, making MB more negative.
    extreme_target = Sleipnir.Float(-5.0)  # -5 m w.e. yr⁻¹ — unreachable by DDF/prcp_fac alone
    dhdt_data = Sleipnir.DhdtData(
        (Sleipnir.Float(calibration_tspan[1]), Sleipnir.Float(calibration_tspan[2])),
        extreme_target)
    glacier_extreme = Sleipnir.Glacier2D(glacier; dhdtData = dhdt_data, geodetic_MB = extreme_target)

    TI_tb = calibrate_ti_model(glacier_extreme, params)
    @test TI_tb isa TImodel1
    # temp_bias must have been used (non-zero) since DDF/prcp_fac couldn't bracket
    @test TI_tb.temp_bias != Sleipnir.Float(0.0)
    # A very negative target requires warming (positive temp_bias)
    @test TI_tb.temp_bias > Sleipnir.Float(0.0)
    mb_tb = compute_mean_annual_MB(
        TI_tb, glacier_extreme,
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(mb_tb)
end

function calibrate_ti_model_default_dhdt_test()
    rgi_ids = ["RGI60-11.03638"]
    calibration_tspan = (2000.0, 2020.0)
    params = _calibration_test_params(rgi_ids, calibration_tspan)

    glacier = Sleipnir.initialize_glaciers(rgi_ids, params)[1]

    @test !isnothing(glacier.dhdtData)
    @test glacier.dhdtData.t == (
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(glacier.dhdtData.dhdt)
    @test glacier.dhdtData.dhdt ≈ Sleipnir.Float(-1.0494)

    TI_cal = calibrate_ti_model(glacier, params)
    @test TI_cal isa TImodel1
    @test TI_cal.DDF >= params.physical.DDF_min
    @test TI_cal.DDF <= params.physical.DDF_max

    mb_cal = compute_mean_annual_MB(
        TI_cal, glacier,
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(mb_cal)
    @test abs(mb_cal - glacier.dhdtData.dhdt) < 1e-4
end
