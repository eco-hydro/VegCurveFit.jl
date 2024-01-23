module VegCurveFit

import Libdl
using Parameters
using Dates
using Statistics

include("tools.jl")

include("CurveFit/CurveFit.jl")
include("optim/optim.jl")
include("weights/wFUN.jl")

end # module
