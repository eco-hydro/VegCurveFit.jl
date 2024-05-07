using VegCurveFit
using Plots
using JLD2, UnPack


@testset "plot_input" begin
  @unpack y, w, t, QC_flag, ylu, dt, nptperyear = load("../data/phenofit_CA-NS6.jld2")
  @test_nowarn p = plot_input(t, y, QC_flag)
end

# write_fig("Figure1--plot_input.pdf", 10, 6)
