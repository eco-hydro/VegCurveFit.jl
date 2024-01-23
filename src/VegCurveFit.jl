module VegCurveFit

using LinearAlgebra
import Libdl
using Parameters
using Dates
using Statistics

include("tools.jl")

include("Optim/Optim.jl")
include("weights/wFUN.jl")

include("get_ylu.jl")

include("CurveFit/CurveFit.jl")

include("smooth_HANT/smooth_HANTS.jl")
include("smooth_SG/main_SG.jl")
include("smooth_whittaker/main_whittaker.jl")

end # module
