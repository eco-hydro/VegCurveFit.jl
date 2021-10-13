function get_bound(lims, keys) 
    lower = [lims[key][1] for key in keys]
    upper = [lims[key][2] for key in keys]
    lower, upper
end

function init_Zhang(par0::param_struct, lims)
    @unpack doy_peak, mn, mx, sos, eos, k, t1, t2 = par0
    prior  = [
        [doy_peak, mn, mx, sos   , k  , eos   , k  ],
        [doy_peak, mn, mx, sos+t1, k*2, eos-t1, k*2],
        [doy_peak, mn, mx, sos-t1, k  , eos+t2, k]]

    keys = ["t0", "mn", "mx", "sos", "k", "eos", "k"]
    prior, get_bound(lims, keys)...
end

function init_AG(par0::param_struct, lims)
    @unpack doy_peak, mn, mx, sos, eos, k, t1, t2, half = par0
    prior = [
        # [doy_peak, mn, mx, 0.2*half, 1  , 0.2*half, 1),
        # [doy_peak, mn, mx, 0.5*half, 1.5, 0.5*half, 1.5),
        [doy_peak, mn, mx, 1/half    , 2, 1/half, 2],
        [doy_peak, mn, mx, 1/(0.8*half), 3, 1/(0.8*half), 3]]
    
    # referenced by TIMESAT
    lower = [lims["t0"][1], lims["mn"][1], lims["mx"][1], 1/(1.4*half), 2, 1/(1.4*half), 2]
    upper = [lims["t0"][2], lims["mn"][2], lims["mx"][2], 1/(0.1*half), 6, 1/(0.1*half), 6]
    prior, lower, upper
end

function init_Beck(par0::param_struct, lims)
    @unpack mn, mx, sos, eos, k, t1, t2= par0
    prior = [
        [mn, mx, sos   , k  , eos   , k], 
        [mn, mx, sos+t1, k*2, eos-t2, k*2]]

    keys = ["mn", "mx", "sos", "k", "eos", "k"]
    prior, get_bound(lims, keys)...
end

function init_Elmore(par0::param_struct, lims)
    @unpack mn, mx, sos, eos, k, t1, t2 = par0
    prior = [
        # [mn, mx - mn, sos   , k*1.25 , eos   , k*1.25 , 0.002),
        # [mn, mx - mn, sos   , k*1    , eos   , k      , 0.001),
        [mn, mx - mn, sos+t1, k*2.5  , eos-t2, k*2.5  , 0.002],
        [mn, mx - mn, sos-t1, k*0.25 , eos+t2, k*0.25 , 0.001]]
    
    # formula = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    # par = ["mn", "mx", "sos", "rsp", "eos", "rau", "m7"]
    keys = ["mn", "mx", "sos", "k", "eos", "k"]
    lower, upper = get_bound(lims, keys)
    prior, [lower; 0], [upper; 1] # add m7
end

function init_Gu(par0::param_struct, lims)
    a  = par0.ampl
    # b1 = 0.1
    # b2 = 0.1
    # c1 = 1
    # c2 = 1
    @unpack mn, mx, sos, eos, k, t1, t2, half = par0
    prior = [
        # [mn, a, a, sos-t1, k/2 , eos+t2, k/2, 1  , 1),
        # [mn, a, a, sos   , k*2 , eos   , k*2, 3  , 3),
        [mn, a, a, sos   , k   , eos   , k  , 2  , 2  ],
        [mn, a, a, sos+t1, k*2 , eos-t2, k*2, 0.5, 0.5],
        [mn, a, a, sos+t1, k*3 , eos-t2, k*3, 5  , 5  ]]
    
    # formula = y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    # par = ["y0", "a1", "a2", "sos", "rsp", "eos", "rau", "c1", "c2"]
    keys = ["mn", "mx", "mx", "sos", "k", "eos", "k"]
    lower, upper = get_bound(lims, keys)
    prior, [lower; 0; 0], [upper; Inf; Inf] # add c1, c2
end

function init_Klos(par0::param_struct, lims)
    a1 = 0
    a2 = 0  #ok
    b1 = par0.mn #ok
    b2 = 0  #ok
    c  = 0.2 * par0.mx  # ok
    ## very slightly smoothed spline to get reliable maximum
    # tmp = smooth.spline(y, df = 0.5 * length(y))

    @unpack doy_peak, mn, mx, sos, eos, k, t1, t2, half = par0
    B1 = 4/(doy_peak - sos)
    B2 = 3.2/(eos - doy_peak)
    m1 = sos + 0.5 * (doy_peak - sos)
    m2 = doy_peak + 0.5 * (eos - doy_peak)
    m1_bis = sos
    m2_bis = eos
    q1 = 0.5  # ok
    q2 = 0.5  # ok
    v1 = 2    # ok
    v2 = 2    # ok  
    prior = [
        [a1, a2, b1, b2, c, B1, B2, m1, m2, q1, q2, v1, v2],
        [a1, a2, b1, 0.01, 0, B1, B2, m1, m2_bis, q1, 1, v1, 4],
        [a1, a2, b1, b2, c, B1, B2, m1_bis, m2, q1, q2, v1, v2],
        [a1, a2, b1, b2, c, B1, B2, m1, m2_bis, q1, q2, v1, v2],
        [a1, a2, b1, b2, c, B1, B2, m1_bis, m2, q1, q2, v1, v2]]
    prior, nothing, nothing
end

function init_Beck(input::input_struct)
    param, lims = init_param(input)
    init_Beck(param, lims)
end

## -----------------------------------------------------------------------------
function FitDL_template(input::input_struct, method = nothing; options...)
    par0, lims = init_param(input)
    
    func_init = getfield(curvefit, Symbol("init_$method"))
    func!     = getfield(curvefit, Symbol("doubleLog_$(method)!"))
    
    prior, lower, upper = func_init(par0, lims)
    # opt (Dict): par, obj
    res = optim_pheno(prior, input, func!;
        lower = lower, upper = upper, options...)
    # Dict("par"    => res["par"], 
    #      "obj"    => res["obj"], 
    #      "method" => method)
end

FitDL_AG(input::input_struct; options...)     = FitDL_template(input, "AG"    ; options...)
FitDL_Zhang(input::input_struct; options...)  = FitDL_template(input, "Zhang" ; options...)
FitDL_Beck(input::input_struct; options...)   = FitDL_template(input, "Beck"  ; options...)
FitDL_Elmore(input::input_struct; options...) = FitDL_template(input, "Elmore"; options...)
FitDL_Gu(input::input_struct; options...)     = FitDL_template(input, "Gu"    ; options...)
FitDL_Klos(input::input_struct; options...)   = FitDL_template(input, "Klos"  ; options...)

function FitDL_Beck2(input::input_struct; options...)
    param, lims = init_param(input)
    prior, lower, upper = init_Beck(param, lims)
    optim_pheno(prior, input, doubleLog_Beck!;
        lower = lower, upper = upper, options...)
end


export FitDL_AG, FitDL_Zhang, FitDL_Beck, FitDL_Elmore, FitDL_Gu, FitDL_Klos, 
    FitDL_Beck2
