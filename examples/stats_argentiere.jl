using Muninn, Sleipnir, Statistics
const RGI = "RGI60-11.03638";
const TSPAN = (2000.0, 2020.0)
rgi_paths = Sleipnir.get_rgi_paths()
params = Sleipnir.Parameters(simulation = Sleipnir.SimulationParameters(
    tspan = TSPAN, multiprocessing = false, workers = 1,
    rgi_paths = Dict(RGI=>rgi_paths[RGI]), use_MB = true, use_velocities = false, climate_data_source = :W5E5))
glacier = Sleipnir.initialize_glaciers([RGI], params)[1]
cal = calibrate_ti_model(glacier, params)
t0, t1 = Sleipnir.Float(TSPAN[1]), Sleipnir.Float(TSPAN[2])
mb = Muninn.compute_cumulative_MB(cal, glacier, t0, t1) ./ (TSPAN[2]-TSPAN[1])
v = mb[.!glacier.mask]
@printf("Muninn (t_melt=0, DDF=%.4f mm/°C/d, n_cells=%d):\n  mean=%.4f  median=%.4f  min=%.4f  max=%.4f m w.e./yr\n",
    cal.DDF*1000, length(v), mean(v), median(v), minimum(v), maximum(v))
