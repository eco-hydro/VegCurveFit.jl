using JLD2
using RData
using UnPack

l = load("data/phenofit-CA_NS6.rda")
@unpack INPUT, brks2 = l

dt = brks2["dt"]
nptperyear = INPUT["nptperyear"] |> Int
y = INPUT["y"]
t = INPUT["t"]
w = INPUT["w"]

ylu = INPUT["ylu"]
ylu = [ylu[1], ylu[2]]

jldsave("data/phenofit-CA_NS6.jld2", true; y, t, w, ylu, nptperyear, dt)
