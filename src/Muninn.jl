__precompile__() # this module is safe to precompile
module Muninn

# ##############################################
# ###########       PACKAGES     ##############
# ##############################################

using Infiltrator
import Pkg
using Distributed
using Dates

### ODINN.jl dependencies ###
using Reexport
using Roots: find_zero, Brent
@reexport using Sleipnir
using Sleipnir: Parameters, Model

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

const src_dir::String = dirname(@__FILE__)
const global root_dir::String = joinpath(src_dir, "..")

# Physical constants for mass balance computation
# Unit conversion: 1 mm of precipitation = 0.001 m w.e.
const PRECIP_UNIT_CONVERSION = Sleipnir.Float(1.0 / 1000.0)
export PRECIP_UNIT_CONVERSION

# All structures and functions related to ODINN models
include(src_dir*"/models/mass_balance/MBmodel.jl")

# Plotting utilities
include(src_dir*"/plotting/plotting_utils.jl")

end # module
