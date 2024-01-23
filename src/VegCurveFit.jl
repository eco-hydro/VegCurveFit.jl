module VegCurveFit

using LinearAlgebra
using Libdl
using Parameters, UnPack
using Dates
using Statistics

using Reexport
@reexport using Serialization: serialize, deserialize
@reexport using DelimitedFiles: readdlm


include("DataTypes.jl")
include("tools.jl")
include("get_ylu.jl")

include("Optim/Optim.jl")
include("weights/wFUN.jl")

include("Smooth/Smooth.jl")
include("CurveFit/CurveFit.jl")
include("VegSeasons/VegSeasons.jl")

end # module
