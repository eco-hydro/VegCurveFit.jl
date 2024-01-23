using Test
using JLD2, UnPack, DataFrames
@unpack y, w, t, ylu, dt, nptperyear = load("../data/phenofit-CA_NS6.jld2")

@testset "curvefits" begin
  @test_nowarn begin
    methods = ["AG", "Zhang", "Beck", "Elmore", "Gu", "Klos"]
    @time r = curvefits(y, t, w, ylu, nptperyear, dt; methods)
  end
end
