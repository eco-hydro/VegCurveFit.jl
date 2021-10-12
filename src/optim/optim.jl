include("nlminb.jl")
include("goal.jl")


"""
    optim_pheno(prior, input; lower = lower, upper = upper, options...)
"""
function optim_pheno(prior, input::input_struct, FUN!::Function = doubleLog_Beck!; 
    lower = nothing, upper = nothing, iters = 2, wFUN = nothing, options...)

    ypred = ones(length(input.y)) .* - 0.1; # initial value of -0.1
    target(par) = goal!(par, FUN!, input, ypred)
    
    ## In each iteration, adjust weights to approach the upper envelope
    opt = nothing;
    for i = 1:iters
        println("i = $i")
        opts = map(par0 -> begin
            nlminb(par0, target; 
                lower = lower, upper = upper, eval_max = 1000, iter_max = 1000)
        end, prior)
        
        # get the best option
        fx, i_opt = map(opt -> opt["objective"], opts) |> findmin
        opt = opts[i_opt]
    end
    opt
end


export optim_pheno;
