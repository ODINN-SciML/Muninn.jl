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
@reexport using Sleipnir

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())
const global root_plots::String = joinpath(root_dir, "plots")

# All structures and functions related to ODINN models
include("models/mass_balance/MBmodel.jl")

end # module

