using Test
using UnPack
using VegCurveFit

@testset "whittaker smoother" begin
  y = [5.0, 8, 9, 10, 12, 10, 15, 10, 9, 19, 19, 17, 13, 14, 18, 19, 18, 12, 18,
    24, 0, 1, 18, 17, 6, 13, 12, 10, 9, 6, 6, 3, 4, 3, 3, 3, 2, 3, 4, 4, 3, 2, 3, 3, 1, 3]
  y = [y; y; y; y]
  m = length(y)
  FT = Float32
  w = ones(Float32, m)

  lambda = 2.0
  z = ones(m)
  interm = interm_whit{FT}(; n=length(y))
  z, cve = whit2!(y, w, lambda, interm; include_cve=true)
  z, cve2 = whit2(y, w; lambda)
  @test cve ≈ cve2

  lamb_cv = lambda_cv(y, w, is_plot=true)
  lamb_vcurve = lambda_vcurve(y, w, is_plot=true)

  z1, cve_cv = whit2(y, w, lambda=lamb_cv)
  z2, cve_vcurve = whit2(y, w, lambda=lamb_vcurve)
  @test cve_cv < cve
  @test cve_vcurve < cve
  @test cve_cv < cve_vcurve
end


@testset "whit2" begin
  d = deserialize("../data/Tumbarumba_EVI2")
  @unpack y, t, w = d

  z1, cve1 = whit2_cv(y, w, lambda=2)
  z2, cve2 = whit2(y, w, lambda=2)
  @test cve1 ≈ cve2
  @test_nowarn r = smooth_whit(y, w)

  # whit3
  x = 1:length(y)
  z1, cve1 = whit3(y, w; lambda=2)
  z2, cve2 = WHIT(y, w, x; lambda=2, d=3)

  @test cve1 ≈ cve2
  @test maximum(z1 - z2) <= 1e-10

  # whit2
  z1, cve1 = whit2(y, w; lambda=2)
  z2, cve2 = WHIT(y, w, x; lambda=2, d=2)

  @test cve1 ≈ cve2
  @test maximum(z1 - z2) <= 1e-10
end

# using BenchmarkTools
# @benchmark
# 4 times faster than R
# @time for i in 1:1e4
#     lamb_vcurve = lambda_vcurve(y, w)
# end
# @time lamb_vcurve = lambda_vcurve(y, w)

# @time @benchmark lamb_cv = lambda_cv(y, w)
# BenchmarkTools.Trial: 
#   memory estimate:  412.55 KiB
#   allocs estimate:  374
#   --------------
#   minimum time:     172.600 μs (0.00% GC)
#   median time:      179.300 μs (0.00% GC)
#   mean time:        198.880 μs (3.88% GC)
#   maximum time:     1.882 ms (70.58% GC)
