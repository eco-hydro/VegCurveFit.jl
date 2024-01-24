function smooth_HANTS(y, w, agrs...;
  FUN = HANTS, wFUN=wBisquare,
  nptperyear=46,
  iters=3,
  param = (nf=3, periodlen=365),
  param_wFUN=(;),
  wmin=0.2, wmid=0.5, wmax=1.0,
  alpha=0.02,
  ignored...)
  
  ylu, wc = get_ylu(y, w; wmin, wmid, wmax, alpha)
  data = DataFrame(; y, w, args...)

  # w1 = ones(size(y));
  res = map(i -> begin
      yfit = FUN(y, w; param...)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu

      w = wFUN(y, yfit, w; iter=i, nptperyear, wmin=0.05, param_wFUN...)
      dfit = DataFrame(; z=yfit, w, agrs...)
    end, 1:iters)

  predict = res # melt_list(res, iter=1:iters)
  Dict("data" => data, "predict" => predict, "param" => param)
end


export HANTS, smooth_HANTS
