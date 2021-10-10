module doubleLogistics
import Libdl
# using Printf

export goal, goal!, gof_RMSE


include("optim/nlminb.jl")
include("f_goal.jl")

include("curveFitting/doubleLog_solve.jl")

end # module
