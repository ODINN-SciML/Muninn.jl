

function TI_creation_default_test(save_refs::Bool = false)

    params = Parameters()
    @inferred Parameters()
    TI1 = TImodel1(params)
    @inferred TImodel1(params)
    TI2 = TImodel2(params)
    @inferred TImodel2(params)

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
    @inferred Parameters()
    TI1 = TImodel1(
        params;
        DDF = 6.0/1000.0,
        acc_factor = 1.2/1000.0
    )
    @inferred TImodel1(
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
    @inferred TImodel2(
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

