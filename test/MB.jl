

function apply_MB_test(save_refs::Bool = false)

    rgi_ids = ["RGI60-11.03638"]

    rgi_paths = get_rgi_paths()
    # Filter out glaciers that are not used to avoid having references that depend on all the glaciers processed in Gungnir
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation = SimulationParameters(
            use_MB=true,
            velocities=false,
            tspan=(2010.0, 2015.0),
            test_mode = true,
            rgi_paths = rgi_paths),
    )
    JET.@test_opt target_modules=(Sleipnir,Muninn) Parameters(
        simulation = SimulationParameters(
            use_MB=true,
            velocities=false,
            tspan=(2010.0, 2015.0),
            test_mode = true,
            rgi_paths = rgi_paths),
    )
    glaciers = initialize_glaciers(rgi_ids, params)
    # JET.@test_opt target_modules=(Sleipnir,Muninn) initialize_glaciers(rgi_ids, params) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    TI1 = TImodel1(params)
    JET.@test_opt TImodel1(params)
    model = Model(nothing, TI1, nothing) # This test only needs a mass balance model
    JET.@test_opt Model(nothing, TI1, nothing)

    glacier = initialize_glaciers(rgi_ids, params)[1]
    t = 2015.0
    step = 1.0/12.0
    mb = MB_timestep(model, glacier, step, t)
    # JET.@test_opt target_modules=(Sleipnir,Muninn) MB_timestep(model, glacier, step, t) # RasterStack manipulation is type unstable, so for the moment this test is deactivated

    if save_refs
        jldsave(joinpath(Muninn.root_dir, "test/data/MB/MB_model.jld2"); mb)
    end

    mb_ref = load(joinpath(Muninn.root_dir, "test/data/MB/MB_model.jld2"))["mb"]
    @test mb == mb_ref

end
