
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

A structure representing a temperature index model with degree-day factor and accumulation factor.

# Keyword arguments
- `DDF::F`: Degree-day factor, which is a coefficient used to convert temperature into melt.
- `acc_factor::F`: Accumulation factor, which is a coefficient used to adjust the accumulation of mass.

# Type Parameters
- `F`: A subtype of `AbstractFloat` representing the type of the factors.
"""
struct TImodel1{F <: AbstractFloat} <: TImodel
    DDF::F
    acc_factor::F
end

"""
    TImodel1(params::Sleipnir.Parameters; DDF::F = 7.0/1000.0, acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

Create a temperature index model with one degree-day factor (DDF) with the given parameters.

# Arguments
- `params::Sleipnir.Parameters`: The simulation parameters.
- `DDF::F`: Degree-day factor (default is 7.0/1000.0).
- `acc_factor::F`: Accumulation factor (default is 1.0/1000.0).

# Returns
- `TI1_model`: An instance of TImodel1 with the specified parameters.
"""
function TImodel1(params::Sleipnir.Parameters;
            DDF::F = 7.0/1000.0,
            acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

    # Build the simulation parameters based on input values
    TI1_model = TImodel1{Sleipnir.Float}(DDF, acc_factor)

    return TI1_model
end

"""
    TImodel2{F <: AbstractFloat}

A type representing a temperature-index model with parameters for snow and ice degree-day factors, and an accumulation factor.

# Keyword arguments
- `DDF_snow::F`: Degree-day factor for snow, which determines the melt rate of snow per degree above the melting point.
- `DDF_ice::F`: Degree-day factor for ice, which determines the melt rate of ice per degree above the melting point.
- `acc_factor::F`: Accumulation factor, which scales the accumulation of snow.

# Type Parameters
- `F`: A subtype of `AbstractFloat`, representing the numeric type used for the model parameters.
"""
struct TImodel2{F <: AbstractFloat} <: TImodel
    DDF_snow::F
    DDF_ice::F
    acc_factor::F
end

"""
    TImodel2(params::Parameters; DDF_snow::F = 3.0/1000.0, DDF_ice::F = 6.0/1000.0, acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

Create a temperature-index model with two degree-day factors (TImodel2) for mass balance calculations.

# Arguments
- `params::Parameters`: The parameters object containing simulation settings.
- `DDF_snow::F`: Degree-day factor for snow (default: 3.0/1000.0).
- `DDF_ice::F`: Degree-day factor for ice (default: 6.0/1000.0).
- `acc_factor::F`: Accumulation factor (default: 1.0/1000.0).

# Returns
- `TI2_model`: An instance of the TImodel2 with the specified parameters.
"""
function TImodel2(params::Parameters;
            DDF_snow::F = 3.0/1000.0,
            DDF_ice::F = 6.0/1000.0,
            acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

    # Build the simulation parameters based on input values
    TI2_model = TImodel2{Sleipnir.Float}(DDF_snow, DDF_ice, acc_factor)

    return TI2_model
end

Base.:(==)(a::TImodel1, b::TImodel1) = a.DDF == b.DDF && a.acc_factor == b.acc_factor

Base.:(==)(a::TImodel2, b::TImodel2) = a.DDF_snow == b.DDF_snow && a.DDF_ice == b.DDF_ice &&
                                        a.acc_factor == b.acc_factor


# Display setup
function Base.show(io::IO, type::MIME"text/plain", model::TImodel1)
    Base.show(io, model)
end
function Base.show(io::IO, model::TImodel1)
    println("Temperature index mass balance model TImodel1")
    print("   DDF = ")
    println(model.DDF)
    print("   acc_factor = ")
    print(model.acc_factor)
end
function Base.show(io::IO, type::MIME"text/plain", model::TImodel2)
    Base.show(io, model)
end
function Base.show(io::IO, model::TImodel2)
    println("Temperature index mass balance model TImodel2")
    print("   DDF_snow = ")
    println(model.DDF_snow)
    print("   DDF_ice = ")
    println(model.DDF_ice)
    print("   acc_factor = ")
    print(model.acc_factor)
end

include("mass_balance_utils.jl")
