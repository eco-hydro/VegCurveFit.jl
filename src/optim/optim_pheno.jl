include("nlminb.jl")
include("goal.jl")


# optim_pheno(prior, sFUN, input; 
#   lower = lower, upper = upper, options...)
function optim_pheno(prior, input::input_struct, FUN!::Function; 
    lower = nothing, upper = nothing, options...)

    ypred = ones(length(input.y)) .* - 0.1; # initial value of -0.1
    # target(par) = goal!(par, FUN!, input, ypred)
    
    par_opts = []
    for i = 1:length(prior)
        par0 = copy(prior[i])
        # @show target(par0)
        par_opt = nlminb(par0, goal!, FUN!, input, ypred; 
            lower = lower, upper = upper, eval_max = 1000, iter_max = 1000)
        # push!(par_opts, par_opt)
    end
    # par_opts
end


export optim_pheno;
