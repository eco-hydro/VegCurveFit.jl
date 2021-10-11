first(x::AbstractArray{T, 1}) where {T<:Real} = x[1] 
last(x::AbstractArray{T, 1}) where {T<:Real} = x[end] 


"""
    init_param(y, t, w = nothing; w_min = 0.5)

Initial parameters for double logistics

# parameters
- `w_min`: weights greater than w_min are treated as good values.

# Examples
```julia
par0, lims = init_param(y, t, w)
```
"""
function init_param(y, t, w = nothing; w_min = 0.5)
    # if (any(is.na(y)))
    #     stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
    if (w === nothing); w = ones(length(y)); end
    
    # t = t - t[1] # (20190103) seq_along(y)  
    # fixed 2018-07-25, If have no enough good points, then set w_min=0
    if (sum(w .>= w_min)/length(y) < .4); w_min = 0; end

    mx     = maximum(y[w .>= w_min])
    mn     = minimum(y[w .>= w_min])
    avg    = mean(y)

    _, pos = findmax(y)
    doy_peak = t[pos]
    # fixed 06 March, 2018; Avoid doy_peak not in the range of doy
    # doy    = quantile(t, c(0.25, 0.75), na.rm = TRUE),
    sos, eos = (doy_peak + first(t))/2, (last(t) + doy_peak) /2
    t1   = (doy_peak - sos)/3 # adjust for sos
    t2   = (eos - doy_peak)/3 # adjust for eos
    # doy  = [sos, eos]
    # if (doy[1] >= doy_peak) doy[1] = (doy_peak - first(t))/2 + first(t)
    # if (doy[2] <= doy_peak) doy[2] = (last(t) - doy_peak) /2 + doy_peak
    ampl   = mx - mn
    deltaY = ampl*0.1
    tmax = maximum(t)
    tmin = minimum(t)
    half   = (tmax - tmin)/2
    deltaT = half/4
    
    k      = 4/half*2.67 #approximate value
    # k limits: about 0.004 - 0.2
    # kmin = 4 / (half * 5), half = 200, k = 0.004
    # kmax = 4 / (half / 5), half = 100, k = 0.2

    # parameters limit
    lims = Dict(
        "t0"  => [doy_peak - deltaT, doy_peak + deltaT],
        "mn"  => [mn - deltaY    , mn + deltaY],
        "mx"  => [mx - deltaY*2  , mx + deltaY*2],
        "r"   => [k/1.2          , k*5],
        "sos" => [tmin         , doy_peak + deltaT],
        "eos" => [doy_peak - deltaT, tmax]
    )
    
    # keys = ["mx", "mn", "ampl", "sos", "eos", "doy_peak", "deltaT", "deltaY", "half", "t1", "t2", "k"]
    # vals = [mx, mn, ampl, sos, eos, doy_peak, deltaT, deltaY, half, t1, t2, k]
    # par0 = Dict(keys .=> vals)
    par0 = param_struct(mx, mn, ampl, sos, eos, doy_peak, deltaT, deltaY, half, t1, t2, k)
    par0, lims
end

function init_param(input::input_struct; options...) 
    init_param(input.y, input.t, input.w; options...)
end

export first, last, init_param
