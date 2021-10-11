include("nlminb.jl")
include("goal.jl")


"""
    optim_pheno(prior, input; lower = lower, upper = upper, options...)
"""
function optim_pheno(prior, input::input_struct, FUN!::Function = doubleLog_Beck!; 
    lower = nothing, upper = nothing, options...)

    ypred = ones(length(input.y)) .* - 0.1; # initial value of -0.1
    # target(par) = goal!(par, FUN!, input, ypred)
    
    par_opts = []
    for i = 1:length(prior)
        par0 = copy(prior[i])
        par_opt = nlminb(par0, goal!, FUN!, input, ypred; 
            lower = lower, upper = upper, eval_max = 1000, iter_max = 1000)
        push!(par_opts, par_opt)
    end
    par_opts
end


export optim_pheno;
