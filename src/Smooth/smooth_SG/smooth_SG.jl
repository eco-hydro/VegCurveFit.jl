function smooth_SG(y, qc, date;
  nptperyear=46,
  iters=3,
  halfwin=5, d=2,
  wmin=0.2, wmid=0.5, wmax=1.0,
  alpha=0.02,
  wFUN=wBisquare,
  options=(;),
  ignored...)

  w, QC_flag = qc_FparLai(qc; wmin=wmin, wmid=wmid, wmax=wmax)
  ylu, wc = get_ylu(y, w; wmin=wmin, wmid=wmid, wmax=wmax, alpha=alpha)

  data = DataFrame(; date, y, w, QC_flag)
  halfwin = round(Int, halfwin)
  # λs = []
  # λᵢ = 2 # default value
  # w1 = ones(size(y));
  res = map(i -> begin
      # if i <= 2
      #   λᵢ = isnothing(λ) ? fun_λ(y, w, is_plot=false) / adj_factor : λ
      # end
      # push!(λs, λᵢ)
      yfit = wSG(y, w; halfwin=halfwin, d=d)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu
      # w[yfit .<= ylu[1]] .= wmin

      w = wFUN(y, yfit, w; iter=i, nptperyear=nptperyear, wmin=0.05, options...)
      dfit = DataFrame(; date, z=yfit, w)
    end, 1:iters)

  predict = melt_list(res, iter=1:iters)
  Dict("data" => data, "predict" => predict,
    "param" => Dict("halfwin" => halfwin, "d" => d))
end


export smooth_SG;
