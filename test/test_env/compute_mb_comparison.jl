using Muninn, Sleipnir, JLD2, Statistics, Printf

# Set up test conditions (same as apply_MB_test)
rgi_ids = ["RGI60-11.03638"]
rgi_paths = Sleipnir.get_rgi_paths()
rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

params = Sleipnir.Parameters(
    simulation = Sleipnir.SimulationParameters(
    use_MB = true,
    use_velocities = false,
    tspan = (2010.0, 2015.0),
    test_mode = true,
    rgi_paths = rgi_paths),
)

@printf("Initializing glacier...\n")
glacier = Sleipnir.initialize_glaciers(rgi_ids, params)[1]

@printf("Creating TI model...\n")
TI1 = Muninn.TImodel1(params)
model = Sleipnir.Model(nothing, TI1, nothing)

@printf("Computing MB timestep...\n")
t = 2015.0
step_MB = 1.0/12.0
mb = Muninn.MB_timestep(model, glacier, step_MB, t)

@printf("\nNew MB statistics:\n")
@printf("  Shape: %s\n", size(mb))
@printf("  Min: %.6f m w.e.\n", minimum(mb))
@printf("  Max: %.6f m w.e.\n", maximum(mb))
@printf("  Mean: %.6f m w.e.\n", mean(mb))
@printf("  Median: %.6f m w.e.\n", median(mb))
@printf("  Std: %.6f m w.e.\n\n", std(mb))

# Load reference
@printf("Loading reference...\n")
mb_ref = load(joinpath(Muninn.root_dir, "test/data/MB/MB_model.jld2"))["mb"]

@printf("Reference MB statistics:\n")
@printf("  Shape: %s\n", size(mb_ref))
@printf("  Min: %.6f m w.e.\n", minimum(mb_ref))
@printf("  Max: %.6f m w.e.\n", maximum(mb_ref))
@printf("  Mean: %.6f m w.e.\n", mean(mb_ref))
@printf("  Median: %.6f m w.e.\n", median(mb_ref))
@printf("  Std: %.6f m w.e.\n\n", std(mb_ref))

# Compute differences
diff = mb - mb_ref
@printf("Difference (NEW - OLD) statistics:\n")
@printf("  Min: %.6f m w.e.\n", minimum(diff))
@printf("  Max: %.6f m w.e.\n", maximum(diff))
@printf("  Mean: %.6f m w.e.\n", mean(diff))
@printf("  Median: %.6f m w.e.\n", median(diff))
@printf("  Std: %.6f m w.e.\n\n", std(diff))

# Save new MB for potential use as new reference
@printf("Saving new MB to temporary file...\n")
jldsave("new_MB_model.jld2"; mb)
@printf("Done.\n")
