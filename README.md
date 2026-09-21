[![Build Status](https://github.com/ODINN-SciML/Muninn.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ODINN-SciML/Muninn.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ODINN-SciML/Muninn.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ODINN-SciML/Muninn.jl)
[![CompatHelper](https://github.com/ODINN-SciML/Muninn.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/ODINN-SciML/Muninn.jl/actions/workflows/CompatHelper.yml)


<img src="https://github.com/JordiBolibar/Muninn.jl/blob/main/data/Muninn_logo-18.png" width="250">

## About Muninn.jl

Muninn.jl is a package providing glacier mass balance models for [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl). It currently includes basic temperature-index models, with 1 or 2 degree-day factors (`TImodel1` and `TImodel2`). In the future, it will also accommodate machine learning models to simulate glacier surface mass balance in a data-driven way, including some physical constraints, see [MassBalanceMachine.jl](https://github.com/ODINN-SciML/MassBalanceMachine.jl).

Temperature-index models can be automatically calibrated per glacier against the geodetic mass balance observations of [Hugonnet et al. (2021)](https://doi.org/10.1038/s41586-021-03436-z) (`calibrate_MB_model`, `compute_mean_annual_MB`), so no parameter needs to be set by hand.

Muninn is part of the ODINN ecosystem, where each package has a narrow role:

  - [Gungnir](https://github.com/ODINN-SciML/Gungnir) (Python): preprocesses OGGM glacier and climate data.
  - [Sleipnir.jl](https://github.com/ODINN-SciML/Sleipnir.jl): core data structures (glaciers, climate, parameters, laws, results).
  - **Muninn**: surface mass balance models (this package).
  - [Huginn.jl](https://github.com/ODINN-SciML/Huginn.jl): ice flow models and PDE solvers.
  - [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl): differentiable pipeline for UDE training and inversions.

## Use Muninn directly or ODINN.jl?

Use Muninn on its own to compute or calibrate mass balance independently of ice flow, or to prototype a new mass balance model. Use [Huginn.jl](https://github.com/ODINN-SciML/Huginn.jl) for coupled ice flow and mass balance forward simulations, and [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl) if you also need gradients through that coupling.

## Installing Muninn

> `Muninn.jl` requires Julia v1.11.

In order to install `Muninn` in a given environment, just do in the REPL:
```julia
julia> ] # enter Pkg mode
(@v1.11) pkg> activate MyEnvironment # or activate whatever path for the Julia environment
(MyEnvironment) pkg> add Muninn
```

Muninn re-exports Sleipnir, so a single `using Muninn` gives access to both packages. The preprocessed glacier data are downloaded automatically the first time Sleipnir is precompiled, see the [Sleipnir README](https://github.com/ODINN-SciML/Sleipnir.jl#data-preprocessing).

## How to use Muninn

The following example calibrates a temperature-index model for two glaciers against the Hugonnet et al. (2021) geodetic mass balance, without setting any parameter by hand:

```julia
using Muninn

# Multiprocessing is disabled for local runs. The Hugonnet et al. observation
# period (2000-2020) is used as the calibration tspan.
params = Parameters(
    simulation = SimulationParameters(
        tspan = (2000.0, 2020.0),
        multiprocessing = false,
        use_MB = true,
        use_velocities = false,
        rgi_paths = get_rgi_paths()
    )
)

# Initializing the glaciers also loads their Hugonnet geodetic mass balance
rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]
glaciers = initialize_glaciers(rgi_ids, params)

# Temperature-index model, with no ice flow component
model = Model(iceflow = nothing, mass_balance = TImodel1(params))

# Calibrate one model per glacier (DDF, prcp_fac and temp_bias are fitted)
model = calibrate_MB_model(model, glaciers, params)

# Compare the observed and the calibrated mean annual mass balance
for (glacier, cal_model) in zip(glaciers, model.mass_balance)
    cal_mb = compute_mean_annual_MB(cal_model, glacier, 2000.0, 2020.0)
    @show glacier.rgi_id, glacier.dhdtData.dhdt, cal_mb # in m w.e. yr⁻¹
end
```

For more details, see the [SMB calibration tutorial](https://odinn-sciml.github.io/ODINN.jl/dev/smb_calibration/). Muninn's own page is [here](https://odinn-sciml.github.io/ODINN.jl/dev/Packages/muninn/), the full list of types and functions is in the [API reference](https://odinn-sciml.github.io/ODINN.jl/dev/API/api_muninn/), and the steps to add a new mass balance model are in [Extending ODINN](https://odinn-sciml.github.io/ODINN.jl/dev/extending/#add-a-new-mass-balance-model).

## Contributing and community

Contributions are welcome. You can report bugs and request features in the [issues](https://github.com/ODINN-SciML/Muninn.jl/issues) tab, or open a pull request against `main` from a fork. See [How to contribute](https://odinn-sciml.github.io/ODINN.jl/dev/contribute/) and the [Code of conduct](https://odinn-sciml.github.io/ODINN.jl/dev/code_of_conduct/) for the guidelines shared across the ODINN ecosystem.

## How to cite

If you use Muninn, please cite the ODINN paper published in [Geoscientific Model Development](https://gmd.copernicus.org/articles/16/6671/2023/gmd-16-6671-2023.html):
```
@article{bolibar_sapienza_universal_2023,
	title = {Universal differential equations for glacier ice flow modelling},
	author = {Bolibar, J. and Sapienza, F. and Maussion, F. and Lguensat, R. and Wouters, B. and P\'erez, F.},
	journal = {Geoscientific Model Development},
	volume = {16},
	year = {2023},
	number = {22},
	pages = {6671--6687},
	url = {https://gmd.copernicus.org/articles/16/6671/2023/},
	doi = {10.5194/gmd-16-6671-2023}
}
```
