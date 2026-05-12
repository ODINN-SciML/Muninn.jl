
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
    glacier = Sleipnir.Glacier2D(glacier; dhdtData = dhdt_data)

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

    # acc_factor should be unchanged (default value)
    @test TI_cal.acc_factor ≈ Sleipnir.Float(1.0 / 1000.0)

    # The modelled MB with the calibrated model should be very close to the target
    mb_cal = compute_mean_annual_MB(
        TI_cal, glacier,
        Sleipnir.Float(calibration_tspan[1]),
        Sleipnir.Float(calibration_tspan[2]))
    @test isfinite(mb_cal)
    @test abs(mb_cal - geodetic_mb) < 1e-4  # within 0.1 mm w.e. yr⁻¹

    # --- Vector method: calibrate over a vector of glaciers ----------------------
    glaciers = [glacier]
    TI_vec = calibrate_ti_model(glaciers, params)
    @test length(TI_vec) == 1
    @test TI_vec[1].DDF ≈ TI_cal.DDF

    # --- calibration_period override --------------------------------------------
    TI_override = calibrate_ti_model(
        glacier, params;
        calibration_period = (Sleipnir.Float(2005.0), Sleipnir.Float(2015.0)))
    @test TI_override isa TImodel1
    @test isfinite(TI_override.DDF)
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
