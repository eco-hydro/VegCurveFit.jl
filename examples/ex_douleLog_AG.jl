# using RCall

# par = R"c( mn  = 0.1, mx  = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)" |> rcopy
# par = Dict(
#     "mn" => 0.1
# )
# t = R"seq(1, 365, 8)"  |> rcopy
# par = R"c(0.1, 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)" |> rget

using doubleLogistics

mn  = 0.1
mx  = 0.7
sos = 50; eos = 250
rsp = 0.1; rau = 0.1

par = [mn, mx, sos, rsp, eos, rau]
t = collect(1.0:8:365)
y = doubleLog_Beck(par, t)
w = ones(length(t))

par0, lims = init_param(y, t, w)


using Plots
plot(t, y)
# solve doubleLog_Beck

struct hello
    x; y
end
