export wSG, whit2, smooth

# include("smooth_HANT/smooth_HANTS.jl")
abstract type AbstractSmoothParam end
include("DataType.jl")

include("HANTS/wHANTS.jl")
include("Whittaker/main_whittaker.jl")
include("SG/main_SG.jl")

include("smooth_whit.jl")
# include("smooth_SG.jl")
# include("smooth_HANTS.jl")


function wSG(y::AbstractVector{FT}, w::AbstractVector{FT}; 
  halfwin=1, d=2, check_wmin=false, kw...) where {FT<:Real}

  z = zeros(FT, size(y))
  wSG!(z, y, w; halfwin, d, check_wmin, kw...) # return z
end

wSG(y::AbstractVector{FT}; kw...) where {FT<:Real} = 
  wSG(y, ones(FT, size(y)); kw...)




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
function smooth(y, w, args...;
  # FUN = wHANTS, params = (nf=3, periodlen=365),
  FUN = wSG, params = (halfwin=5, d=2),
  nptperyear=46, iters=3,
  wmin::Float32=0.2f0, wmid::Float32=0.5f0, wmax::Float32=1.0f0,
  
  alpha=0.02,
  wFUN=wBisquare,
  options=(;),
  # use_spike=false,
  ignored...)
  
  # w, QC_flag = qc_FparLai(qc; wmin, wmid, wmax)
  ylu, wc = get_ylu(y, w; wmin, wmid, wmax, alpha)
  data = DataFrame(; y, w, args...)

  # λs = []
  # λᵢ = 2 # default value
  res = map(i -> begin
      # auto adjust parameters
      # push!(λs, λᵢ)
      yfit = FUN(y, w; params...)
      # yfit, cve = whit2(y, w, λᵢ)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu
      # w[yfit .<= ylu[1]] .= wmin
      w = wFUN(y, yfit, w; iter=i, nptperyear, wmin=0.05, options...)
      dfit = DataFrame(; date, y0=y, z=yfit, w)
    end, 1:iters)

  predict = res # melt_list(res, iter=1:iters)
  # (; data, predict, param=λs)
  Dict("data" => data, "predict" => predict, param => params)
end
