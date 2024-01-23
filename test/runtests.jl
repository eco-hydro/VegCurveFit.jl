using Test
using VegCurveFit

@testset "nlminb" begin
  # Write your own tests here.
  start = [0.0] # [1, 2, 3, 4]
  lower = [-4.0]
  upper = [4.0]

  function f2(x::Array{Float64,1})::Float64
    x[1] - cos(x[1])
  end

  f(x) = x[1] - cos(x[1])
  r = nlminb(start, f, verbose=false)
  @test r["par"][1] == -1.5701006351106073
end
