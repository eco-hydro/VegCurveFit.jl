# include("smooth_HANT/smooth_HANTS.jl")
include("smooth_HANTS.jl")
include("smooth_SG/main_SG.jl")
include("smooth_whittaker/main_whittaker.jl")


"""
    smooth(y, qc, date;
        niters=5,
        is_plot=true, title="whittaker",
        outfile="Figures/Plot-smooth_whit.pdf",
        options...)
    
# Arguments
- `alpha`     : for quantile calculating in `get_ylu`
- `options...`: other parameters to [wBisquare()]
"""
function smooth(y, w, QC_flag, date;
  iters=3,
  wmin::Float32=0.2f0, wmid::Float32=0.5f0, wmax::Float32=1.0f0,
  nptperyear=46,
  alpha=0.02,
  wFUN=wBisquare,
  options=(;),
  # use_spike=false,
  ignored...)

  # w, QC_flag = qc_FparLai(qc; wmin, wmid, wmax)
  ylu, wc = get_ylu(y, w; wmin, wmid, wmax, alpha)
  data = DataFrame(; date, y, w, QC_flag)

  # λs = []
  # λᵢ = 2 # default value
  res = map(i -> begin
      # auto adjust parameters
      # push!(λs, λᵢ)
      yfit = FUN(y, w)
      # yfit, cve = whit2(y, w, λᵢ)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu
      # w[yfit .<= ylu[1]] .= wmin
      w = wFUN(y, yfit, w; iter=i, nptperyear, wmin=0.05, options...)
      dfit = DataFrame(; date, y0=y, z=yfit, w)
    end, 1:iters)

  predict = res # melt_list(res, iter=1:iters)
  # (; data, predict, param=λs)
  Dict("data" => data, "predict" => predict)
end
