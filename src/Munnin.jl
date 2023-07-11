__precompile__() # this module is safe to precompile
module Munnin

# ##############################################
# ###########       PACKAGES     ##############
# ##############################################

using Infiltrator
import Pkg
using Distributed
using PyCall
using Sleipnir

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())
const global root_plots::String = joinpath(root_dir, "plots")

# All structures and functions related to ODINN models
include(joinpath(ODINN.root_dir, "src/models/MBmodel.jl"))

end # module

