module curvefit

import Libdl
using Parameters
using Dates
using Statistics

include("tools.jl")

include("curveFitting/curvefits.jl")
include("optim/optim.jl")
include("weights/wFUN.jl")

end # module
