
export compute_MB, MB_timestep!, MB_timestep,
       requires_dynamic_topography, topography_window_m, mb_inputs,
       required_climate_data_source, validate_climate_data_source,
       validate_model_simulation_compatibility

requires_dynamic_topography(::MBmodel) = false
topography_window_m(::MBmodel) = Sleipnir.Float(200.0)
mb_inputs(::MBmodel) = (;)
required_climate_data_source(::MBmodel) = nothing

function validate_climate_data_source(model::MBmodel, climate_data_source::Symbol)
    required_source = required_climate_data_source(model)
    if !isnothing(required_source) && climate_data_source != required_source
        throw(ArgumentError(
            "Mass-balance model $(typeof(model)) requires climate_data_source = :$(required_source), got :$(climate_data_source)."
        ))
    end
    return nothing
end

function validate_model_simulation_compatibility(
        model::Sleipnir.Model,
        parameters::Sleipnir.Parameters)
    if !isnothing(model.mass_balance)
        validate_climate_data_source(
            model.mass_balance,
            parameters.simulation.climate_data_source
        )
    end
    return nothing
end

function _requires_topography_from_inputs(model::MBmodel)
    inputs = values(mb_inputs(model))
    return any(inp -> inp isa Union{Sleipnir.iTopoSlope, Sleipnir.iTopoAspect}, inputs)
end

"""
    compute_MB(
        mb_model::TImodel1,
        climate_2D_period::Climate2Dstep,
        step::AbstractFloat,
    )

Compute the mass balance (MB) for a given mass balance model and climate period.

# Arguments

  - `mb_model::TImodel1`: The mass balance model containing parameters such as accumulation factor (`acc_factor`) and degree-day factor (`DDF`).
  - `climate_2D_period::Climate2Dstep`: The climate data for a specific period, including snow accumulation (`snow`) and positive degree days (`PDD`).
  - `step::AbstractFloat`: The step used to update MB. This scales the MB so that the accumulation and degree-day factors are scaled monthly.

# Returns

  - A numerical array representing the computed mass balance, calculated as the difference between the product of the accumulation factor and snow, and the product of the degree-day factor and positive degree days.
"""
function compute_MB(
        mb_model::TImodel1,
        climate_2D_period::Climate2Dstep,
        step::AbstractFloat
)
    return ((mb_model.acc_factor .* mb_model.prcp_fac .* climate_2D_period.snow) .-
            (mb_model.DDF .* climate_2D_period.PDD)) ./ (step / (1/12))
end

"""
    MB_timestep(model::Model, glacier::G, step::F, t::F) where {F <: AbstractFloat, G <: AbstractGlacier}

Calculate the mass balance (MB) for a glacier over a given timestep.

# Keyword arguments

  - `model::Model`: The model containing mass balance parameters.
  - `glacier::G`: The glacier object containing climate and DEM data.
  - `step::F`: The timestep duration.
  - `t::F`: The current time.

# Returns

  - `MB::Matrix{F}`: The computed mass balance matrix for the given timestep.

# Details

 1. Computes the period between the current time `t` and the previous step `t - step`.
 2. Retrieves cumulative climate data for the specified period.
 3. Downscales the climate data to a 2D grid based on the glacier's DEM.
 4. Computes the mass balance using the downscaled climate data.
"""
function MB_timestep(model::Model, glacier::G, step::F,
        t::F) where {F <: AbstractFloat, G <: AbstractGlacier}
    get_cumulative_climate!(glacier.climate, t, step)

    # Convert climate dataset to 2D based on the glacier's DEM
    climate_2D_step = downscale_2D_climate(
        glacier.climate.climate_step,
        glacier.S,
        glacier.Coords;
        include_topography = requires_dynamic_topography(model.mass_balance) ||
                             _requires_topography_from_inputs(model.mass_balance),
        topography_window_m = topography_window_m(model.mass_balance),
        Δx = glacier.Δx,
        Δy = glacier.Δy)
    return compute_MB(model.mass_balance, climate_2D_step, step)
end

"""
    MB_timestep!(cache, model::Model, glacier::G, step::F, t) where {F <: AbstractFloat, G <: AbstractGlacier}

Simulates a mass balance timestep for a given glacier model.

# Arguments

  - `cache`: The model cache to update.
  - `model::Model`: The glacier model containing iceflow and mass balance information.
  - `glacier::G`: The glacier object containing climate and DEM data.
  - `step::F`: The timestep duration.
  - `t`: The current time.

# Description

This function performs the following steps:

 1. Computes the period from the previous timestep to the current time.
 2. Retrieves cumulative climate data for the specified period.
 3. Downscales the climate dataset to a 2D grid based on the glacier's DEM.
 4. Computes the mass balance for the glacier and updates the model's iceflow mass balance.
"""
function MB_timestep!(cache, model::Model, glacier::G, step::F,
        t) where {F <: AbstractFloat, G <: AbstractGlacier}
    get_cumulative_climate!(glacier.climate, t, step)

    # Convert climate dataset to 2D based on the glacier's DEM
    downscale_2D_climate!(
        glacier;
        include_topography = requires_dynamic_topography(model.mass_balance) ||
                             _requires_topography_from_inputs(model.mass_balance),
        topography_window_m = topography_window_m(model.mass_balance))
    cache.iceflow.MB .= compute_MB(
        model.mass_balance, glacier.climate.climate_2D_step, step)
    return nothing # For type stability
end
