import Statistics: quantile


get_range(y::AbstractVector{T}) where {T<:Real} = [minimum(y), maximum(y)]

get_ylu(y) = get_range(y)


"""
    get_ylu(y::AbstractVector{<:Real}, w::AbstractVector{<:Real};
        wmin=0.2, wmid=0.5, wmax=1.0, 
        ratio = 0.4, 
        alpha=0.02, alpha_high=nothing)
"""
function get_ylu(y::AbstractVector{<:Real}, w::AbstractVector{<:Real};
  wmin=0.2, wmid=0.5, wmax=1.0,
  ratio=0.4,
  alpha=0.02, alpha_high=nothing)

  alpha_high = null_default(alpha_high, alpha)

  w_critical = wmin
  n = length(y)

  if (sum(w .== wmax) >= n * ratio)
    w_critical = wmax
  elseif (sum(w .>= wmid) > n * ratio)
    # Just set a small portion for boreal regions. In this way, it will give
    # more weights to marginal data.
    w_critical = wmid
  end
  y_good = @view y[w.>=w_critical]

  # alpha/2, alpha_high set to 0.05 for remote sensing (20200322)
  ymin = max(quantile(y_good, alpha / 2), 0) # TODO: fix CN-Din
  ymax = quantile(y_good, 1 - alpha_high / 2)
  ylu = [ymin, ymax]
  ylu, w_critical
end


export get_ylu;
