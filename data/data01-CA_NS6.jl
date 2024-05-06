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


# QC_flag = INPUT["QC_flag"]
n = length(y)
w = ones(n)
t = 1:n

# 快20倍左右
z, cve = whit3(y, w; lambda=2);
# cve_right = 0.03960781157141451
@time z, cve = whit3(y, w; lambda=2);
@time z2, cve2 = WHIT(y, t, w; d=3, lambda=2);
@test maximum(z - z2) <= 1e-10
@test abs(cve - cve2) <= 1e-10

## whit2
@time z, cve = whit2(y, w; lambda=2)
@time z2, cve2 = WHIT(y, t, w; d=2, lambda=2)
cve - cve2

@btime z, cve, h = whit3($y, $w; lambda=2);
@btime z2 = WHIT($y, $t; d=3, lambda=2);

@btime z, cve = whit2($y, $w; lambda=2);
@btime z2 = WHIT($y, $t; d=2, lambda=2);

begin
  # z, cve = whit2(y, ones(n); lambda=2)
  # z2 = whit(y, 1:length(y); d=2, lambda=2)
  # @time z, cve = whit3(y, ones(n); lambda=2)
  # @time z2 = whit(y, 1:length(y); d=3, lambda=2)
  # z - z2
  using Plots
  plot()
  # plot(y, label="y", color="grey")
  plot!(z, label="fast")
  plot!(z2, label="matrix")
end

# ylu = INPUT["ylu"]
# ylu = [ylu[1], ylu[2]]
# using RTableTools
# d = DT(; y, w)
# fwrite(d, "data.csv")
# jldsave("data/phenofit-CA_NS6.jld2", true; y, t, w, QC_flag, ylu, nptperyear, dt)
