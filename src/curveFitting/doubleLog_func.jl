using Statistics


"""
doubleLog_AG
"""
function doubleLog_AG(par::T, t::T)::T where T <: AbstractArray{<:Real, 1}
    # if (sos >= eos) return(rep(9999, length(t)))    
    ypred = ones(T, length(t)) * 99.0
    doubleLog_AG!(ypred, par, t); ypred
end

function doubleLog_AG!(ypred::T, par::T, t::T) where T <: AbstractArray{<:Real, 1}
    t0  = par[1]
    mn  = par[2]
    mx  = par[3]
    rsp = par[4]
    a3  = par[5]
    rau = par[6]
    a5  = par[7]

    # if (sos >= eos) return(rep(9999, length(t)))    
    n = length(t)
    @inbounds for i = 1:n
        if (t[i] <= t0)
            ypred[i] = mn + (mx - mn)*exp(- ((t0 - t[i])*rsp) ^a3 )
        else
            ypred[i] = mn + (mx - mn)*exp(  ((t[i] - t0)*rau) ^a5 ) 
        end
    end
    # ypred = mn + (mx - mn)*exp(- c( ((t0 - t[t <= t0])*rsp) ^a3,
    #                 ((t[t >  t0] - t0)*rau) ^a5) )
end


"""
doubleLog_Beck
"""
function doubleLog_Beck(par::AbstractArray{T, 1}, t::AbstractArray{T, 1}) where T <: Real
    ypred = ones(T, length(t)) * 99.0
    doubleLog_Beck!(ypred, par, t); ypred
end

function doubleLog_Beck!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T, 1}) where T <: Real
    # mn  = par[1]
    # mx  = par[2]
    sos = par[3]
    eos = par[5]
    rsp = par[4]
    rau = par[6]
    
    if (eos < sos)
        ypred .= 99.0
    else
        for i = 1:length(t)
            ypred[i] = par[1] + (par[2] - par[1]) * (
                1 / (1 + exp(-rsp * (t[i] - sos))) + 
                1 / (1 + exp(rau * (t[i] - eos))) - 1 )
        end
    end
end

export doubleLog_AG, doubleLog_AG!, doubleLog_Beck, doubleLog_Beck!
