using Test
using JLD2, UnPack
@unpack y, w, t, ylu, dt, nptperyear = load("data/phenofit-CA_NS6.jld2")

@testset "curvefits" begin
  @test_nowarn begin
    @time r = curvefits(y, t, w, ylu, nptperyear, dt)
  end
end
