
function TI_creation_default_test(save_refs::Bool = false)
    params = Parameters()
    JET.@test_opt target_modules=(Sleipnir, Muninn) Parameters()
    TI1 = TImodel1(params)
    JET.@test_opt TImodel1(params)
    TI2 = TImodel2(params)
    JET.@test_opt TImodel2(params)

    if save_refs
        jldsave(joinpath(Muninn.root_dir, "test/data/TI/TI1_model_default.jld2"); TI1)
        jldsave(joinpath(Muninn.root_dir, "test/data/TI/TI2_model_default.jld2"); TI2)
    end

    TI1_ref = load(joinpath(Muninn.root_dir, "test/data/TI/TI1_model_default.jld2"))["TI1"]
    TI2_ref = load(joinpath(Muninn.root_dir, "test/data/TI/TI2_model_default.jld2"))["TI2"]

    @test TI1 == TI1_ref
    @test TI2 == TI2_ref
end

function TI_creation_values_test(save_refs::Bool = false)
    params = Parameters()
    JET.@test_opt target_modules=(Sleipnir, Muninn) Parameters()
    TI1 = TImodel1(
        params;
        DDF = 6.0/1000.0,
        acc_factor = 1.2/1000.0
    )
    JET.@test_opt TImodel1(
        params;
        DDF = 6.0/1000.0,
        acc_factor = 1.2/1000.0
    )
    TI2 = TImodel2(
        params;
        DDF_snow = 3.0/1000.0,
        DDF_ice = 6.0/1000.0,
        acc_factor = 1.2/1000.0
    )
    JET.@test_opt TImodel2(
        params;
        DDF_snow = 3.0/1000.0,
        DDF_ice = 6.0/1000.0,
        acc_factor = 1.2/1000.0
    )

    if save_refs
        jldsave(joinpath(Muninn.root_dir, "test/data/TI/TI1_model_specified.jld2"); TI1)
        jldsave(joinpath(Muninn.root_dir, "test/data/TI/TI2_model_specified.jld2"); TI2)
    end

    TI1_ref = load(joinpath(Muninn.root_dir, "test/data/TI/TI1_model_specified.jld2"))["TI1"]
    TI2_ref = load(joinpath(Muninn.root_dir, "test/data/TI/TI2_model_specified.jld2"))["TI2"]

    @test TI1 == TI1_ref
    @test TI2 == TI2_ref
end

function TI_synthetic_field_test()
    params = Parameters()
    F = Sleipnir.Float

    snow = F[2.0 0.5; 0.0 1.0]
    PDD = F[0.0 1.0; 3.0 0.25]
    z2 = zeros(F, 2, 2)
    xy = F[0.0, 1.0]

    climate_2D_step = Sleipnir.Climate2Dstep(
        temp = z2,
        PDD = PDD,
        snow = snow,
        rain = z2,
        elevation_diff = z2,
        aspect = z2,
        albedo = z2,
        slhf = z2,
        slope = z2,
        sshf = z2,
        ssrd = z2,
        str = z2,
        gradient = F(0.0),
        avg_gradient = F(0.0),
        x = xy,
        y = xy,
        ref_hgt = F(0.0)
    )

    TI1 = TImodel1(
        params;
        DDF = F(4.0/1000.0),
        acc_factor = F(2.0/1000.0)
    )

    # Explicit TI equation terms to avoid opaque hardcoded references.
    accumulation_term = TI1.acc_factor .* snow
    melt_term = TI1.DDF .* PDD

    step_month = F(1.0/12.0)
    scale_month = step_month / F(1.0/12.0)
    expected_month = (accumulation_term .- melt_term) ./ scale_month

    mb_month = compute_MB(TI1, climate_2D_step, step_month)
    @test mb_month == expected_month
    @test expected_month == F[0.004 -0.003; -0.012 0.001]

    # The implementation rescales MB by step/(1/12), so halving the step doubles MB.
    step_half = F(1.0/24.0)
    scale_half = step_half / F(1.0/12.0)
    expected_half = (accumulation_term .- melt_term) ./ scale_half
    mb_half = compute_MB(TI1, climate_2D_step, step_half)
    @test mb_half == expected_half
    @test expected_half == 2 .* expected_month
end
