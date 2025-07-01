
export compute_MB, MB_timestep!, MB_timestep


"""
    compute_MB(mb_model::TImodel1, climate_2D_period::Climate2Dstep)

Compute the mass balance (MB) for a given mass balance model and climate period.

# Arguments
- `mb_model::TImodel1`: The mass balance model containing parameters such as accumulation factor (`acc_factor`) and degree-day factor (`DDF`).
- `climate_2D_period::Climate2Dstep`: The climate data for a specific period, including snow accumulation (`snow`) and positive degree days (`PDD`).

# Returns
- A numerical array representing the computed mass balance, calculated as the difference between the product of the accumulation factor and snow, and the product of the degree-day factor and positive degree days.
"""
function compute_MB(mb_model::TImodel1, climate_2D_period::Climate2Dstep)
    return ((mb_model.acc_factor .* climate_2D_period.snow) .- (mb_model.DDF .* climate_2D_period.PDD))
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
function MB_timestep(model::Model, glacier::G, step::F, t::F) where {F <: AbstractFloat, G <: AbstractGlacier}
    # First we get the dates of the current time and the previous step
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)
    get_cumulative_climate!(glacier.climate, period)
    # Convert climate dataset to 2D based on the glacier's DEM
    climate_2D_step::Climate2Dstep = downscale_2D_climate(glacier.climate.climate_step, glacier.S, glacier.Coords)
    MB::Matrix{F} = compute_MB(model.mass_balance, climate_2D_step)
    return MB
end


"""
    MB_timestep!(model::Model, glacier::G, step::F, t; batch_id::Union{Nothing, I} = nothing) where {I <: Integer, F <: AbstractFloat, G <: AbstractGlacier}

Simulates a mass balance timestep for a given glacier model.

# Arguments
- `cache`: The model cache to update.
- `model::Model`: The glacier model containing iceflow and mass balance information.
- `glacier::G`: The glacier object containing climate and DEM data.
- `step::F`: The timestep duration.
- `t`: The current time.
- `batch_id::Union{Nothing, I}`: Optional batch identifier for simulations using Reverse Diff. Defaults to `nothing`.

# Description
This function performs the following steps:
1. Computes the period from the previous timestep to the current time.
2. Retrieves cumulative climate data for the specified period.
3. Downscales the climate dataset to a 2D grid based on the glacier's DEM.
4. Computes the mass balance for the glacier and updates the model's iceflow mass balance.

If `batch_id` is provided, the function updates the mass balance for the specified batch; otherwise, it updates the mass balance for the entire model.
"""
function MB_timestep!(cache, model::Model, glacier::G, step::F, t; batch_id::Union{Nothing, I} = nothing) where {I <: Integer, F <: AbstractFloat, G <: AbstractGlacier}
    # First we get the dates of the current time and the previous step
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)

    get_cumulative_climate!(glacier.climate, period)

    # Convert climate dataset to 2D based on the glacier's DEM
    downscale_2D_climate!(glacier)
    # Simulations using Reverse Diff require an iceflow and mass balance model per glacier
    if isnothing(batch_id)
        cache.iceflow.MB .= compute_MB(model.mass_balance, glacier.climate.climate_2D_step)
    else
        cache.iceflow[batch_id].MB .= compute_MB(model.mass_balance[batch_id], glacier.climate.climate_2D_step)
    end
    return nothing # For type stability
end
