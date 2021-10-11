include("data_types.jl")
include("init_param.jl")
include("doubleLog_func.jl")


function FitDL_Beck(input; options...)
    # if (missing(w)) w = rep(1, length(y))
    # input = input_struct(y, t, w)
    par0, lims = init_param(input)
    # sFUN  = "doubleLog_Beck!"
    @unpack mn, mx, sos, eos, k, t1, t2, k = par0
    prior = [
        [mn, mx, sos   , k  , eos   , k], 
        [mn, mx, sos+t1, k*2, eos-t2, k*2]]
    keys = ["mn", "mx", "sos", "r", "eos", "r"]
    lower = [lims[key][1] for key in keys]
    upper = [lims[key][2] for key in keys]
    
    optim_pheno(prior, input, doubleLog_Beck!;
        lower = lower, upper = upper, options...)
end


export FitDL_Beck
