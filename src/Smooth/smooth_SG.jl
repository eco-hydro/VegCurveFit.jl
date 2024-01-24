function smooth_SG(y, w, args...;
  nptperyear=46,
  iters=3,
  halfwin=5, d=2,
  wmin=0.2, wmid=0.5, wmax=1.0,
  alpha=0.02,
  wFUN=wBisquare,
  options=(;),
  ignored...)

  ylu, wc = get_ylu(y, w; wmin, wmid, wmax, alpha)

  data = DataFrame(; y, w, args...)
  halfwin = round(Int, halfwin)
  # w1 = ones(size(y));
  res = map(i -> begin
      yfit = wSG(y, w; halfwin=halfwin, d=d)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu
      # w[yfit .<= ylu[1]] .= wmin

      w = wFUN(y, yfit, w; iter=i, nptperyear=nptperyear, wmin=0.05, options...)
      dfit = DataFrame(; z=yfit, w, args...)
    end, 1:iters)

  param = Dict("halfwin" => halfwin, "d" => d)
  predict = res # melt_list(res, iter=1:iters)
  Dict("data" => data, "predict" => predict, "param" => param)
end


export smooth_SG;
