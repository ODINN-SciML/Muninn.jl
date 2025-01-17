
export compute_MB, MB_timestep!, MB_timestep

function compute_MB(mb_model::TImodel1, climate_2D_period::Climate2Dstep)
    return ((mb_model.acc_factor .* climate_2D_period.snow) .- (mb_model.DDF .* climate_2D_period.PDD))
end

function MB_timestep(model::Model, glacier::G, step::F, t::F) where {F <: AbstractFloat, G <: AbstractGlacier}
    # First we get the dates of the current time and the previous step
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)
    climate_step::Dict{String, Any} = get_cumulative_climate(glacier.climate[At(period)])
    # Convert climate dataset to 2D based on the glacier's DEM
    climate_2D_step::Climate2Dstep = downscale_2D_climate(climate_step, glacier)
    MB::Matrix{F} = compute_MB(model.mb_model, climate_2D_step)
    return MB
end


function MB_timestep!(model::Model, glacier::G, step::F, t; batch_id::Union{Nothing, I} = nothing) where {I <: Integer, F <: AbstractFloat, G <: AbstractGlacier}
    # First we get the dates of the current time and the previous step
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)

    get_cumulative_climate!(glacier.climate, period)

    # Convert climate dataset to 2D based on the glacier's DEM
    downscale_2D_climate!(glacier)
    # Simulations using Reverse Diff require an iceflow and mass balance model per glacier
    if isnothing(batch_id)
        model.iceflow.MB .= compute_MB(model.mass_balance, glacier.climate.climate_2D_step)
    else
        model.iceflow[batch_id].MB .= compute_MB(model.mass_balance[batch_id], glacier.climate.climate_2D_step)
    end
end

