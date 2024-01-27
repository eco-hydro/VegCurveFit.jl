using Dates

@testset "wHANTS" begin
  d = deserialize("../data/Tumbarumba_EVI2")
  @unpack y, t, w = d
  x = datetime2julian.(DateTime.(t)) |> x -> x .- x[1]

  @test_nowarn z2, amp, phi = wHANTS(y, w, x; nf=3)
end
