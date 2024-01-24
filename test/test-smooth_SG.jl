using Test

@testset "Savitzky Golay filter" begin
  # works
  @test wSG([2., 3, 5, 7]; halfwin=1, d=2) ≈ [2, 3, 5, 7]
  @test wSG_low([2., 3, 5, 7]; halfwin=1, d=2) ≈ [2, 3, 5, 7]

  y = 1.0:10 |> collect
  @test wSG(y; halfwin=1, d=1) ≈ y

  # 跟R语言版本基本一致
  y = [2.0, 3, 4, 10, 6, 7]
  w = [1.0, 1, 1, 0.2, 1, 1]
  z1 = wSG(y, w; halfwin=1, d=1)
  z2 = wSG_low(y, w; halfwin=1, d=1)
  
  @test round.(z1, digits=2) ≈ [2, 3, 5, 5.45, 7, 6.5]
  @test round.(z2, digits=2) ≈ [2, 3, 5, 5.45, 7, 6.5]

  # wSG对权重调整没有这么敏感
  z = wSG(y, w; halfwin=1, d=2)
  @test round.(z, digits=2) ≈ y
end


@testset "Savitzky Golay filter" begin
  y = [1.0, 2, 5, 4, 3, 6]
  w = collect(1.0:7)
  halfwin = 3
  S = sgmat_S(halfwin)
  B = sgmat_B(S)
  sgmat_wB(S, w)
  z = SG(y; halfwin=1)
  @test z ≈ y
end


@testset "Savitzky Golay filter" begin
  y = rand(100)
  w = rand(100)
  SG(y; halfwin=5)
  wSG(y, w; halfwin=5)
end


# using Plots
# p = plot(y)
# plot!(p, z1)
# plot!(p, z2)
# n = Int(1e5)
# y = rand(n)
# using BenchmarkTools

## test the used memory
# @time for i=1:1e3
#     z = SG(y)
# end
# using Plots
# gr()
# plot(y)
# plot!(z)
