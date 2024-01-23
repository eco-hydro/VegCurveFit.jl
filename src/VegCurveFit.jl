module VegCurveFit

using LinearAlgebra
using Libdl
using Parameters
using Dates
using Statistics

include("DataTypes.jl")
include("tools.jl")
include("get_ylu.jl")

include("Optim/Optim.jl")
include("weights/wFUN.jl")

include("Smooth/Smooth.jl")
include("CurveFit/CurveFit.jl")
include("VegSeasons/VegSeasons.jl")

end # module
