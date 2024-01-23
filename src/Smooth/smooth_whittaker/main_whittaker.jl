using Statistics
import Parameters: @with_kw, @with_kw_noshow


# include("smooth_whittaker/whit2_cpp.jl")
# include("smooth_whittaker/whittaker2.jl")
# include("smooth_whittaker/lambda_cv.jl")
# include("smooth_whittaker/smooth_whit_v01.jl")
include("lambda_init.jl")
include("lambda_vcurve.jl")
include("whit2.jl")
include("whit2_Frasso2015.jl")

include("smooth_whit.jl")
# include("smooth_whit_GEE.jl")
