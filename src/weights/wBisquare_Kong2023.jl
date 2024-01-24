"""
    wBisquare(y, yfit, w; iter = 2, wmin, to_upper = true)

# Arguments
- `iter`    : not used
- `options` : currently ignored

# Bad points
- 1. under the yfit, in the growing season (yfit > 0.3 * A + ymin)

# Examples
"""
function wBisquare_Kong2023(y::AbstractVector{<:Real}, yfit::AbstractVector{<:Real}, w::AbstractVector{T2};
  # QC_flag;
  # iter::Integer = 2, 
  # to_upper = true, 
  step::Float64=0.3,
  wmin::Float64=0.05, wmax=2.0,
  trs_high=0.7,
  trs_low=0.4,
  trs_bg=0.2,
  verbose=false,
  ignored...) where {T2<:Real}

  # change wBisquare to consider yearly amplitude differences
  ymin = minimum(yfit)
  ymax = maximum(yfit)
  A = ymax - ymin

  re = yfit .- y
  re_abs = abs.(re)
  # println("re: $re")
  sc = 6 * median(re_abs)
  # println("sc", median(re_abs), ", ", sc)
  # 最保险的方法，获取每年的ylu，然后判断是ingrowing or ungrowing
  # println("threshold：high=", trs_high * A + ymin, ", low=", trs_low * A + ymin)
  y_high = trs_high * A + ymin
  y_low = trs_low * A + ymin
  y_bg = trs_bg * A + ymin
  verbose && println("y_high = $y_high, y_low = $y_low, y_bg = $y_bg")

  I_high = @.(yfit > y_high)
  I_low = @.((yfit < y_low) & (yfit > y_bg))

  ## method 2: 强力加强季节性信号
  # y_2 = [0; 0; diff(diff(yfit))]
  # I_high = y_2 .< 0
  # I_low = y_2 .> 0
  ## 这一步出现了问题
  I_bad_high = @.((re > 0) & I_high & (w < 1)) # middle of GS, upper envelope
  I_bad_low = @.((re < 0) & I_low & (w < 1))

  # `I_good_high`可以考虑增加`w > 0.2`的限制
  I_good_high = @.((re < 0) & I_high)
  I_good_low = @.((re > 0) & I_low & (w > 0.2))

  I_bad = I_bad_high .| I_bad_low
  I_good = I_good_high .| I_good_low

  wnew = deepcopy(w)
  wnew[I_bad] = @.((1 - (re_abs[I_bad] / sc)^2)^2 * w[I_bad]) # - step    
  wnew[I_good] = wnew[I_good] .+ step

  ## 3 异常值赋予最低权重，改写对应的value
  # outlier不可纵容，即使是re < 0
  # 相信QC，质量好的点不惩罚

  sc_ext = min(A / 2, sc * 2)
  I_outlier = @.(((re_abs >= sc) & (w <= 0.5) & I_low) |
                 ((re_abs >= sc_ext) & (w < 1)) | 
                 (re_abs >= A))

  wnew[I_outlier] .= wmin
  wnew[wnew.<wmin] .= wmin
  wnew[wnew.>wmax] .= wmax
  
  ## fix `44_CH-Dav`
  y[I_bad_high] .= yfit[I_bad_high] # y was changed at here
  y[I_outlier] .= ymin  
  ## also need to change y to improve the upper envelope performance
  # println("inside wBisquare: ", sum(y))
  # wmax = 2.0
  # wnew[wnew .>= 1.0] .= 1.0
  wnew
end
