# using RCall

# par = R"c( mn  = 0.1, mx  = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)" |> rcopy
# par = Dict(
#     "mn" => 0.1
# )
# t = R"seq(1, 365, 8)"  |> rcopy
# par = R"c(0.1, 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)" |> rget
using Pkg
# Pkg.activate(".")
Pkg.activate("/mnt/e/Research/julia/nlminb.jl")

using doubleLogistics
using Parameters

begin
    mn  = 0.1
    mx  = 0.7
    sos = 50; eos = 250
    rsp = 0.1; rau = 0.1

    par = [mn, mx, sos, rsp, eos, rau]
    t = collect(1.0:8:365)
    y = doubleLog_Beck(par, t)
    w = ones(length(t))
    
    input = input_struct(y, t, w)
    # FitDL_Beck(input)    

    par0, lims = init_param(input)
    # sFUN  = "doubleLog_Beck!"
    @unpack mn, mx, sos, eos, k, t1, t2, k = par0
    prior = [
        [mn, mx, sos   , k  , eos   , k], 
        [mn, mx, sos+t1, k*2, eos-t2, k*2]]
    keys = ["mn", "mx", "sos", "r", "eos", "r"]
    lower = [lims[key][1] for key in keys]
    upper = [lims[key][2] for key in keys]

    ypred = ones(length(y)) * -1
    par0 = [0.1, 0.6, 40, 0.05,200, 0.12]
end

f(x) = x[1] - cos(x[1]) 

# println(nlminb([1.0], f))
@show nlminb([ 1.0 ], f; lower = [-10], upper = [12], verbose = true)



# using Plots
# plot(t, y)
# for i = 1:1
    # ans = goal!(par0, doubleLog_Beck!, input, ypred)
# ans = nlminb(par0, goal!, doubleLog_Beck!, input, ypred);
# print(ans);

# ans = nlminb(par0, goal!, doubleLog_Beck!, input, ypred; 
#     lower = lower, upper = upper)
# print(ans);

#     println("ans = $ans")
# # end




# ## ex1
# sumsq(x, y) = sum((x -y).^2)
# y = repeat([1], 5)
# x0 = rand(5)
# nlminb(x0, sumsq, y)
