module nlminb

import Libdl
# using Printf

export optim_nlminb, goal, goal!, gof_RMSE, 
    doubleLog_Beck, doubleLog_Beck!

include("optim_nlminb.jl")
include("doubleLogistics.jl")
include("f_goal.jl")

end # module
