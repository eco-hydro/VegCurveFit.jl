using Parameters

include("doubleLog_func.jl")
include("init_param.jl")
include("data_types.jl")


function FitDL_Beck(y, t, tout = t, w = nothing; options...)
    # if (missing(w)) w = rep(1, length(y))
    input = input_struct(y, t, w, tout)
    par0, lims = init_param(y, t, w)

    sFUN  = "doubleLog.Beck"
    
    @unpack mn, mx, sos, eos, k, t1, t2, k = par0
    prior = [
        [mn, mx, doy[1]   , k  , doy[2]   , k], 
        [mn, mx, doy[1]+t1, k*2, doy[2]-t2, k*2]]
    keys = ["mn", "mx", "sos", "r", "eos", "r"]
    lower = [lims[key][1] for key in keys]
    upper = [lims[key][1] for key in keys]
    # optim_pheno(prior, sFUN, input; 
    #     lower = lower, upper = upper, options...)
end
