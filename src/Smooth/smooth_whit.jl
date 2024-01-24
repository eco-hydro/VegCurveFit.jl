"""
    smooth_whit(y, qc, date;
        niters=5,
        λ=nothing,
        fun_λ=lambda_cv,
        adj_factor=1.0,
        is_plot=true, title="whittaker",
        outfile="Figures/Plot-smooth_whit.pdf",
        options...)
    
# Arguments

- `adj_factor`: λ = lambda_opt / adj_factor
- `alpha`     : for quantile calculating in `get_ylu`
- `options...`: other parameters to [wBisquare()]
"""
function smooth_whit(y, w, args...;
  iters=3,
  λ=nothing,
  fun_λ=lambda_cv,
  adj_factor=1.0,
  wmin::Float32=0.2f0, wmid::Float32=0.5f0, wmax::Float32=1.0f0,
  nptperyear=46,
  alpha=0.02,
  wFUN=wBisquare,
  options=(;),
  # use_spike=false,
  ignored...)

  # w, QC_flag = qc_FparLai(qc; wmin, wmid, wmax)
  ylu, wc = get_ylu(y, w; wmin, wmid, wmax, alpha)
  data = DataFrame(; y, w, args)

  λs = []
  λᵢ = 2 # default value
  res = map(i -> begin
      if i <= 2
        λᵢ = isnothing(λ) ? fun_λ(y, w, is_plot=false) / adj_factor : λ
      end
      push!(λs, λᵢ)

      yfit, cve = whit2(y, w, λᵢ)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu
      # w[yfit .<= ylu[1]] .= wmin
      w = wFUN(y, yfit, w; iter=i, nptperyear, wmin=0.05, options...)
      dfit = DataFrame(; y0=y, z=yfit, w, args...)
    end, 1:iters)

  predict = res # melt_list(res, iter=1:iters)
  # (; data, predict, param=λs)
  Dict("data" => data, "predict" => predict, "param" => λs)
end

# if use_spike
#   ## not tested, kdd
#   # 2023-08-08, Yuxuan Xie
#   I_spike, ymov = check_spike(y; nptperyear)
#   I_bad = @. (w != wmax) # just change point without `good` target
#   I_bad_spike = @.((I_spike) & (I_bad))
#   y[I_bad_spike] .= minimum(y)
# end

export smooth_whit
