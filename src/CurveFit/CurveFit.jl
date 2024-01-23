include("_get_ylu.jl")
include("init_param.jl")
include("doubleLogistics.jl")

include("get_extent.jl")
include("FitDL.jl")
include("Constant.jl")


"""
# Arguments
- `y::AbstractVector{T}`: input time series

"""
function curvefits(y::AbstractVector{T}, t, w::AbstractVector{T},
  ylu::Vector{Float64}, nptperyear::Integer,
  dt;
  methods=["AG", "Zhang", "Beck", "Elmore"]) where {T<:Real}

  # ylu  = input_R[:ylu]
  input = input_struct(y, t, w)
  # input = input_struct(input_R[:y], input_R[:t], input_R[:w])
  t = input.t
  w = input.w
  w0 = copy(w)
  doys = (t - DATE_ORIGIN) |> day2num
  years = Dates.year.(t)
  # nptperyear = input_R[:nptperyear] |> Int
  width_ylu = nptperyear * 2

  begs = getDateId(dt[:, :beg], t)
  # peaks = getDateId(dt[:, :peak], t)
  ends = getDateId(dt[:, :_end], t, "after")
  ns = length(begs) # number of seasons

  opts = []
  for i = 1:ns
    I_beg = begs[i]
    I_end = ends[i]
    I = I_beg:I_end
    I_extend = get_extentI(w0, I_beg, I_end, nptperyear)
    ylu_period = _get_ylu(input.y, years, w0, width_ylu, I; wmin=0.2) # ylu_period 
    merge_ylu!(ylu, ylu_period)
    # println(ylu_period)

    ti = doys[I_extend]
    yi = input.y[I_extend]
    wi = input.w[I_extend]
    # @show I_extend, ti, yi, wi, ylu_period

    ibeg = i == 1 ? 1 : 2
    tout = ti[ibeg]:ti[end]

    ## solve double logistics
    inputI = input_struct(yi, ti, wi)
    opt = Dict{String,Any}()
    
    for meth = methods
      temp = FUNCS_FITDL[meth](inputI) # par, obj, method
      # opt["ypred"] = getfield(curvefit, Symbol("doubleLog_$meth"))(opt["par"], tout)
      temp["zs"] = FUNCS_doubleLog[meth](temp["par"], tout)
      temp["fun"] = "doubleLog.$meth"
      opt[meth] = temp
      # push!(opt, temp)
    end
    # pars = map(x -> x["par"], opt)
    # pars = pmap(opt, "par")
    # predictor = predictInput_struct(models, tout, ylu)
    # predictor = predictor_struct(opt, tout, ylu)
    r = Dict("tout" => tout, "model" => opt)
    push!(opts, r)
  end
  opts
end


export curvefits
