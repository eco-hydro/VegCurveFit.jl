using Test
using VegCurveFit

@testset "findpeaks" begin
  d = findpeaks([1, 2, 3, 2, 0])
  @test size(d, 1) == 1
  @test d.pos_peak[1] == 3

  d = findpeaks([1, 2, 3, 3, 2, 0])
  @test size(d, 1) == 1
  @test d.pos_peak[1] == 3
end
