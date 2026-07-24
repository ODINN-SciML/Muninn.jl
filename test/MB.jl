mutable struct fakeIceflowCache{F <: AbstractFloat}
    MB::Matrix{F}
end

mutable struct fakeCache{ICEFLOW}
    iceflow::ICEFLOW
end

struct DummyMBDefault <: Muninn.MBmodel end
struct DummyMBEra5 <: Muninn.MBmodel end
struct DummyMBPlainInputs <: Muninn.MBmodel end
struct DummyMBTopoInputs <: Muninn.MBmodel end

Muninn.required_climate_data_source(::DummyMBEra5) = :ERA5
Muninn.mb_inputs(::DummyMBPlainInputs) = (;)
Muninn.mb_inputs(::DummyMBTopoInputs) = (; slope = Sleipnir.iTopoSlope())

function model_compatibility_utils_test()
    params_era5 = Parameters(
        simulation = Sleipnir.SimulationParameters(
        climate_data_source = :ERA5, multiprocessing = false))
    params_w5e5 = Parameters(
        simulation = Sleipnir.SimulationParameters(
        climate_data_source = :W5E5, multiprocessing = false))

    # Default MB model imposes no climate source constraint.
    @test isnothing(validate_climate_data_source(DummyMBDefault(), :W5E5))

    # MB model with explicit climate source requirement validates both paths.
    @test isnothing(validate_climate_data_source(DummyMBEra5(), :ERA5))
    @test_throws ArgumentError validate_climate_data_source(DummyMBEra5(), :W5E5)

    # Compatibility check should be a no-op when there is no MB model.
    model_without_mb = Model(nothing, nothing, nothing)
    @test isnothing(validate_model_simulation_compatibility(model_without_mb, params_w5e5))

    # Compatibility check should enforce required climate source when MB exists.
    model_era5_mb = Model(nothing, DummyMBEra5(), nothing)
    @test isnothing(validate_model_simulation_compatibility(model_era5_mb, params_era5))
    @test_throws ArgumentError validate_model_simulation_compatibility(model_era5_mb, params_w5e5)

    # Topography input detection branch coverage.
    @test !Muninn._requires_topography_from_inputs(DummyMBPlainInputs())
    @test Muninn._requires_topography_from_inputs(DummyMBTopoInputs())
end

function apply_MB_test(save_refs::Bool = false)
    rgi_ids = ["RGI60-11.03638"]

    rgi_paths = get_rgi_paths()
    # Filter out glaciers that are not used to avoid having references that depend on all the glaciers processed in Gungnir
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation = SimulationParameters(
        use_MB = true,
        use_velocities = false,
        tspan = (2010.0, 2015.0),
        test_mode = true,
        multiprocessing = false,
        rgi_paths = rgi_paths),
    )
    JET.@test_opt target_modules=(Sleipnir, Muninn) Parameters(
        simulation = SimulationParameters(
        use_MB = true,
        use_velocities = false,
        tspan = (2010.0, 2015.0),
        test_mode = true,
        multiprocessing = false,
        rgi_paths = rgi_paths),
    )
    glacier = initialize_glaciers(rgi_ids, params)[1]
    # JET.@test_opt broken=true target_modules=(Sleipnir,Muninn) initialize_glaciers(rgi_ids, params) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    TI1 = TImodel1(params)
    JET.@test_opt TImodel1(params)
    model = Model(nothing, TI1, nothing) # This test only needs a mass balance model
    JET.@test_opt Model(nothing, TI1, nothing)

    t = 2015.0
    step_MB = 1.0/12.0
    mb = MB_timestep(model, glacier, step_MB, t)
    JET.@test_opt target_modules=(Sleipnir, Muninn) MB_timestep(model, glacier, step_MB, t)

    iceflowCache = fakeIceflowCache{Sleipnir.Float}(zero(glacier.H₀))
    cache = fakeCache{typeof(iceflowCache)}(iceflowCache)

    MB_timestep!(cache, model, glacier, step_MB, t)
    @assert mb==cache.iceflow.MB
    JET.@test_opt target_modules=(Sleipnir, Muninn) MB_timestep!(
        cache, model, glacier, step_MB, t)

    if save_refs
        jldsave(joinpath(Muninn.root_dir, "test/data/MB/MB_model.jld2"); mb)
    end

    mb_ref = load(joinpath(Muninn.root_dir, "test/data/MB/MB_model.jld2"))["mb"]
    @test mb == mb_ref
end
