include("DataTypes.jl")
include("init_param.jl")
include("doubleLog_func.jl")

include("get_ylu.jl")
include("get_extent.jl")
include("FitDL.jl")


function curvefits(input_R, brks)
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
    peaks = getDateId(brks[:dt][:, :peak], t)
    ends  = getDateId(brks[:dt][:, :end], t, "after")
    ns = length(begs) # number of seasons

    opts = []
    @time for i = 1:ns
        println("i = $i")
        I_beg = begs[i]
        I_end = ends[i]
        I = I_beg:I_end
        I_extend = get_extentI(w0, I_beg, I_end, nptperyear)
        ylu = get_ylu(input.y, years, w0, width_ylu, I; wmin = 0.2)

        ti = doys[I_extend]
        yi = input.y[I_extend]
        wi = input.w[I_extend]
        
        # ibeg = i == 1 ? 1 : 2
        # tout = ti[ibeg]:ti[end]
        ## solve double logistics
        inputI = input_struct(yi, ti, wi)
        opt = FitDL_Beck(inputI)
        push!(opts, opt)
    end
    opts
end


export curvefits
