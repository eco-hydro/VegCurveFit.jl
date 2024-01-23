using nlminb
using Test

@testset "nlminb.jl" begin
  # Write your own tests here.
  start = [0.0] # [1, 2, 3, 4]
  lower = [-4.0]
  upper = [4.0]

  function f2(x::Array{Float64,1})::Float64
    x[1] - cos(x[1])
  end

  f(x) = x[1] - cos(x[1])
  r = optim_nlminb(start, f, verbose=false)
end
