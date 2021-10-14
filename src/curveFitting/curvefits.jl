include("DataTypes.jl")
include("init_param.jl")
include("doubleLogistics.jl")

include("get_ylu.jl")
include("get_extent.jl")
include("FitDL.jl")


FINE_fitings = Dict(
    "AG"     => FitDL_AG,
    "Zhang"  => FitDL_Zhang, 
    "Beck"   => FitDL_Beck,
    "Elmore" => FitDL_Elmore, 
    "Gu"     => FitDL_Gu,
    "Klos"   => FitDL_Klos
)

function curvefits(input_R, brks; 
    methods = ["AG", "Zhang", "Beck", "Elmore"]) # , "Gu"

    ylu  = input_R[:ylu]
    input = input_struct(input_R[:y], input_R[:t], input_R[:w])
    date_origin = Date("2000-01-01")
    t    = input.t
    w    = input.w
    w0   = copy(w)
    doys = (t - date_origin) |> day2num
    years = Dates.year.(t)

    nptperyear = input_R[:nptperyear] |> Int
    width_ylu  = nptperyear*2

    begs  = getDateId(brks[:dt][:, :beg], t)
    # peaks = getDateId(brks[:dt][:, :peak], t)
    ends  = getDateId(brks[:dt][:, :end], t, "after")
    ns = length(begs) # number of seasons

    opts = []
    for i = 1:ns
        I_beg = begs[i]
        I_end = ends[i]
        I = I_beg:I_end
        I_extend = get_extentI(w0, I_beg, I_end, nptperyear)
        ylu_period = get_ylu(input.y, years, w0, width_ylu, I; wmin = 0.2) # ylu_period 
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
        opt = map(meth -> FINE_fitings[meth](inputI), methods)
        # par, obj, method
        
        # pars = map(x -> x["par"], opt)
        # pars = pmap(opt, "par")
        # predictor = predictInput_struct(models, tout, ylu)
        predictor = predictor_struct(opt, tout, ylu)
        push!(opts, predictor)
        # push!(opts, opt)
    end
    opts
end



export curvefits
