
export TImodel1, TImodel2, MBmodel

"""
    MBmodel <: AbstractModel

An abstract type representing a mass balance model in the Muninn package.
This serves as a base type for all specific mass balance models, ensuring
they adhere to a common interface and can be used interchangeably within
the ODINN framework.
"""
abstract type MBmodel <: AbstractModel end

###############################################
########## TEMPERATURE-INDEX MODELS ###########
###############################################

"""
    TImodel <: MBmodel

An abstract type representing a temperature index mass balance models within the ODINN framework.
This type serves as a parent type for more specialized mass balance models, ensuring they adhere to
a common interface defined by the `MBmodel` abstract type.
"""
abstract type TImodel <: MBmodel end

"""
    TImodel1{F <: AbstractFloat}

A structure representing a temperature index model with degree-day factor
and precipitation correction factor.

# Keyword arguments

  - `DDF::F`: Degree-day factor (m w.e. °C⁻¹ d⁻¹), which converts positive degree days to melt.
  - `prcp_fac::F`: Dimensionless precipitation correction factor applied as a multiplier to
    the snowfall field before computing accumulation.  A value of `1.0` (default) leaves
    precipitation unchanged; values > 1 increase accumulation (useful when the climate input
    underestimates solid precipitation), values < 1 reduce it.

# Type Parameters

  - `temp_bias::F`: Uniform temperature bias (°C) added to the glacier's climate
    before computing melt and snow/rain partitioning.  `0.0` (default) leaves the
    climate unchanged.  Used as the third calibration lever when `DDF` and `prcp_fac`
    alone cannot bracket the geodetic mass-balance target.

  - `F`: A subtype of `AbstractFloat` representing the type of the factors.

Note: The unit conversion from mm to m w.e. is handled internally via the constant
`PRECIP_UNIT_CONVERSION` (1/1000) and is not a tunable parameter.
"""
struct TImodel1{F <: AbstractFloat} <: TImodel
    DDF::F
    prcp_fac::F
    temp_bias::F
end

"""
    TImodel1(params::Sleipnir.Parameters; DDF::F = 7.0/1000.0, prcp_fac::F = 1.0, temp_bias::F = 0.0) where {F <: AbstractFloat}

Create a temperature index model with one degree-day factor (DDF) with the given parameters.

# Arguments

  - `params::Sleipnir.Parameters`: The simulation parameters.
  - `DDF::F`: Degree-day factor in m w.e. °C⁻¹ d⁻¹ (default is 7.0/1000.0 = 0.007).
  - `prcp_fac::F`: Dimensionless precipitation correction factor (default is 1.0).
  - `temp_bias::F`: Uniform temperature bias in °C (default is 0.0).

# Returns

  - `TI1_model`: An instance of TImodel1 with the specified parameters.

Note: Precipitation unit conversion (mm → m w.e.) is handled internally via `PRECIP_UNIT_CONVERSION`.
"""
function TImodel1(::Sleipnir.Parameters;
        DDF::F = 7.0/1000.0,
        prcp_fac::F = 1.0,
        temp_bias::F = 0.0) where {F <: AbstractFloat}
    return TImodel1{Sleipnir.Float}(DDF, prcp_fac, temp_bias)
end

"""
    TImodel2{F <: AbstractFloat}

A type representing a temperature-index model with parameters for snow and ice degree-day factors.

# Keyword arguments

  - `DDF_snow::F`: Degree-day factor for snow in m w.e. °C⁻¹ d⁻¹, which determines the melt rate of snow per degree above the melting point.
  - `DDF_ice::F`: Degree-day factor for ice in m w.e. °C⁻¹ d⁻¹, which determines the melt rate of ice per degree above the melting point.

# Type Parameters

  - `F`: A subtype of `AbstractFloat`, representing the numeric type used for the model parameters.

Note: The unit conversion from mm to m w.e. is handled internally via the constant
`PRECIP_UNIT_CONVERSION` (1/1000) and is not a tunable parameter.
"""
struct TImodel2{F <: AbstractFloat} <: TImodel
    DDF_snow::F
    DDF_ice::F
end

"""
    TImodel2(params::Sleipnir.Parameters; DDF_snow::F = 3.0/1000.0, DDF_ice::F = 6.0/1000.0) where {F <: AbstractFloat}

Create a temperature-index model with two degree-day factors (TImodel2) for mass balance calculations.

# Arguments

  - `params::Sleipnir.Parameters`: The parameters object containing simulation settings.
  - `DDF_snow::F`: Degree-day factor for snow in m w.e. °C⁻¹ d⁻¹ (default: 3.0/1000.0 = 0.003).
  - `DDF_ice::F`: Degree-day factor for ice in m w.e. °C⁻¹ d⁻¹ (default: 6.0/1000.0 = 0.006).

# Returns

  - `TI2_model`: An instance of the TImodel2 with the specified parameters.

Note: Precipitation unit conversion (mm → m w.e.) is handled internally via `PRECIP_UNIT_CONVERSION`.
"""
function TImodel2(params::Sleipnir.Parameters;
        DDF_snow::F = 3.0/1000.0,
        DDF_ice::F = 6.0/1000.0) where {F <: AbstractFloat}

    # Build the simulation parameters based on input values
    TI2_model = TImodel2{Sleipnir.Float}(DDF_snow, DDF_ice)

    return TI2_model
end

function Base.:(==)(a::TImodel1, b::TImodel1)
    a.DDF == b.DDF && a.prcp_fac == b.prcp_fac && a.temp_bias == b.temp_bias
end

function Base.:(==)(a::TImodel2, b::TImodel2)
    a.DDF_snow == b.DDF_snow && a.DDF_ice == b.DDF_ice
end

# Display setup
Base.show(io::IO, type::MIME"text/plain", model::TImodel1) = Base.show(io, model)
function Base.show(io::IO, model::TImodel1)
    println(io, "Temperature index mass balance model TImodel1")
    print(io, "   DDF = ")
    println(io, model.DDF)
    print(io, "   prcp_fac = ")
    println(io, model.prcp_fac)
    print(io, "   temp_bias = ")
    print(io, model.temp_bias)
end
Base.show(io::IO, type::MIME"text/plain", model::TImodel2) = Base.show(io, model)
function Base.show(io::IO, model::TImodel2)
    println(io, "Temperature index mass balance model TImodel2")
    print(io, "   DDF_snow = ")
    println(io, model.DDF_snow)
    print(io, "   DDF_ice = ")
    println(io, model.DDF_ice)
end

include("mass_balance_utils.jl")
include("calibration.jl")
