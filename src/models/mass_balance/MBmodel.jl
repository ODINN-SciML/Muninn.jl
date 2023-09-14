
export TImodel1, TImodel2, MBmodel

# Abstract type as a parent type for Mass Balance models
abstract type MBmodel <: AbstractModel end 

###############################################
########## TEMPERATURE-INDEX MODELS ###########
###############################################

# Subtype structure for Temperature-Index Mass Balance model
abstract type TImodel <: MBmodel end
# Temperature-index model with 1 melt factor
# Make these mutable if necessary
struct TImodel1{F <: AbstractFloat} <: TImodel
    DDF::F
    acc_factor::F
end

"""
    TImodel1(params::Parameters;
        DDF::Float64 = 5.0/1000.0,
        acc_factor::Float64 = 1.0/1000.0
        )
Temperature-index model with a single degree-day factor.

Keyword arguments
=================
    - `DDF`: Single degree-day factor, for both snow and ice.
    - `acc_factor`: Accumulation factor
"""
function TImodel1(params::Parameters;
            DDF::F = 5.0/1000.0,
            acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

    # Build the simulation parameters based on input values
    ft = params.simulation.float_type
    TI1_model = TImodel1{ft}(DDF, acc_factor)

    return TI1_model
end

# Temperature-index model with 2 melt factors
struct TImodel2{F <: AbstractFloat} <: TImodel
    DDF_snow::F
    DDF_ice::F
    acc_factor::F
end

"""
    TImodel2(params::Parameters;
        DDF_snow::Float64 = 3.0/1000.0,
        DDF_ice::Float64 = 6.0/1000.0,
        acc_factor::Float64 = 1.0/1000.0
        )
Temperature-index model with two melt factors, for snow and ice.

Keyword arguments
=================
    - `DDF_snow`: Degree-day factor for snow.
    - `DDF_ice`: Degree-day factor for ice.
    - `acc_factor`: Accumulation factor
"""
function TImodel2(params::Parameters;
            DDF_snow::F = 3.0/1000.0,
            DDF_ice::F = 6.0/1000.0,
            acc_factor::F = 1.0/1000.0) where {F <: AbstractFloat}

    # Build the simulation parameters based on input values
    ft = params.simulation.float_type
    TI2_model = TImodel2{ft}(DDF_snow, DDF_ice, acc_factor)

    return TI2_model
end

Base.:(==)(a::TImodel1, b::TImodel1) = a.DDF == b.DDF && a.acc_factor == b.acc_factor 

Base.:(==)(a::TImodel2, b::TImodel2) = a.DDF_snow == b.DDF_snow && a.DDF_ice == b.DDF_ice && 
                                        a.acc_factor == b.acc_factor 

include("mass_balance_utils.jl")
