module VegCurveFit

using LinearAlgebra
using Libdl
using Parameters, UnPack
using Dates
using Statistics

using Reexport
@reexport using Serialization: serialize, deserialize
@reexport using DelimitedFiles: readdlm
@reexport using UnPack

# using Plots
using RecipesBase
using RecipesBase: plot, plot!
@shorthands scatter
@shorthands vline
function stroke end

include("DataTypes.jl")
include("tools.jl")
include("get_ylu.jl")

include("plot_input.jl")

include("Optim/Optim.jl")
include("weights/weights.jl")

include("Smooth/Smooth.jl")
include("CurveFit/CurveFit.jl")
include("VegSeasons/VegSeasons.jl")

end # module
