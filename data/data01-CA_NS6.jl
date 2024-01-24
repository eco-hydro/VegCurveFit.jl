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
QC_flag = INPUT["QC_flag"]

ylu = INPUT["ylu"]
ylu = [ylu[1], ylu[2]]

jldsave("data/phenofit-CA_NS6.jld2", true; y, t, w, QC_flag, ylu, nptperyear, dt)
