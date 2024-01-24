using Test
using VegCurveFit
using DataFrames, Plots

include("test-weights.jl")
include("test-plot_input.jl")

include("test-VegSeasons.jl")
include("test-curvefit.jl")
include("test-lambda_init.jl")
include("test-smooth_SG.jl")
include("test-smooth_whit.jl")
include("test-Optim.jl")
