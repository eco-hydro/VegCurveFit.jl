using JLD2
using RData
using UnPack
using VegCurveFit
using BenchmarkTools
using Test

l = load("data/phenofit-CA_NS6.rda")
@unpack INPUT, brks2 = l

dt = brks2["dt"]
nptperyear = INPUT["nptperyear"] |> Int
y = INPUT["y"]
w = INPUT["w"]
t = INPUT["t"]

serialize("data/phenofit_CA-NS6", (;y, w, t))

# QC_flag = INPUT["QC_flag"]
n = length(y)
w = ones(n)
t = 1:n

# ylu = INPUT["ylu"]
# ylu = [ylu[1], ylu[2]]
# using RTableTools
# d = DT(; y, w)
# fwrite(d, "data.csv")
# jldsave("data/phenofit-CA_NS6.jld2", true; y, t, w, QC_flag, ylu, nptperyear, dt)
