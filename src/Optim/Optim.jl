include("nlminb.jl")
include("goal.jl")


"""
    optim_pheno(prior, input; lower = lower, upper = upper, options...)

# Parameters
- `options`: Other parameters to wFUN

"""
function optim_pheno(prior, input::input_struct, FUN!::Function=doubleLog_Beck!;
  lower=nothing, upper=nothing, iters=2, (wFUN!)=wTSM!, options...)

  ypred = ones(length(input.y)) .* -0.1 # initial value of -0.1
  # target(par) = goal!(par, FUN!, input, ypred)
  w = copy(input.w)
  ## In each iteration, adjust weights to approach the upper envelope
  opt = nothing
  for i = 1:iters
    # println("i = $i")
    opts = map(par0 -> begin
        nlminb(par0, goal!, FUN!, input, ypred;
          lower, upper, feval_max=1000, iter_max=1000)
      end, prior)
    # get the best option
    fx, i_opt = map(opt -> opt["obj"], opts) |> findmin
    opt = opts[i_opt]
    # update weight in each iteration
    FUN!(ypred, opt["par"], input.t)
    wFUN!(input.y, ypred, input.w; iter=i, nptperyear=23, wfact=0.5) # update weight by reference
  end
  input.w = w
  opt
end


export optim_pheno;
