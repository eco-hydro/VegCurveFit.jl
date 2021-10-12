"""
doubleLog_AG
"""
function doubleLog_AG!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    t0  = par[1]
    mn  = par[2]
    mx  = par[3]
    rsp = par[4]
    a3  = par[5]
    rau = par[6]
    a5  = par[7]

    # if (sos >= eos) return(rep(9999, length(t)))    
    @inbounds for i = 1:length(t)
        if (t[i] <= t0)
            ypred[i] = mn + (mx - mn)*exp(- ((t0 - t[i])*rsp) ^a3 )
        else
            ypred[i] = mn + (mx - mn)*exp(  ((t[i] - t0)*rau) ^a5 ) 
        end
    end
end

function doubleLog_Zhang!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    t0  = par[1]
    mn  = par[2]
    mx  = par[3]
    sos = par[4]
    rsp = par[5]
    eos = par[6]
    rau = par[7]

    if (t0 - sos <= 1 || t0 - eos >= -1);  ypred .= 99.0; end

    @inbounds for i = 1:length(t)
        if (t[i] <= t0)
            ypred[i] = mn + (mx - mn)/(1 + exp(-rsp*(t[i] - sos)))
        else
            ypred[i] = mn + (mx - mn)/(1 + exp( rau*(t[i] - eos)))
        end
    end
end

"""
doubleLog_Beck
"""
function doubleLog_Beck!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    # mn  = par[1]
    # mx  = par[2]
    sos = par[3]
    eos = par[5]
    rsp = par[4]
    rau = par[6]
    
    for i = 1:length(t)
        ypred[i] = par[1] + (par[2] - par[1]) * (
            1 / (1 + exp(-rsp * (t[i] - sos))) + 
            1 / (1 + exp(rau * (t[i] - eos))) - 1 )
    end
end

function doubleLog_Elmore!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    # mn  = par[1]
    # mx  = par[2]
    sos = par[3] # SOS
    rsp = par[4] # 1/rsp
    eos = par[5] # EOS
    rau = par[6] # 1/rau
    m7  = par[7]
    for i = 1:length(t)
        ypred[i] = par[1] + (par[2] - m7*t[i])*( 
            1/(1 + exp(-rsp*(t[i]-sos))) - 
            1/(1 + exp(-rau*(t[i]-eos))) ) 
    end
end

function doubleLog_Gu!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    y0  = par[1]
    a1  = par[2]
    a2  = par[3]
    sos = par[4]
    rsp = par[5]
    eos = par[6]
    rau = par[7]
    c1  = par[8]
    c2  = par[9]
    for i = 1:length(t)
        # ypred = y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
        ypred[i] = y0 + (a1/(1 + exp(-rsp*(t[i] - sos)))^c1) - (a2/(1 + exp(-rau*(t[i] - eos)))^c2)
    end
end


function doubleLog_Klos!(ypred::AbstractArray{T, 1}, par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real}
    a1 = par[1]
    a2 = par[2]
    b1 = par[3]
    b2 = par[4]
    c  = par[5]
    B1 = par[6]
    B2 = par[7]
    m1 = par[8]
    m2 = par[9]
    q1 = par[10]
    q2 = par[11]
    v1 = par[12]
    v2 = par[13]
    for i = 1:length(t)
        ypred[i] = (a1*t[i] + b1) + (a2*t[i]^2 + b2*t[i] + c) * 
            (1/(1 + q1 * exp(-B1 * (t[i] - m1)))^v1 - 1/(1 + q2 * exp(-B2 * (t[i] - m2)))^v2)
    end
end

## FUNCTIONS for USERS ---------------------------------------------------------

function doubleLog_template(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}, FUN!) where {T<:Real, T2<:Real}
    ypred = ones(T, length(t)) * 99.0
    FUN!(ypred, par, t)
    ypred
end

doubleLog_AG(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_AG!)

doubleLog_Zhang(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_Zhang!)

doubleLog_Beck(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_Beck!)

doubleLog_Elmore(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_Elmore!)

doubleLog_Gu(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_Gu!)

doubleLog_Klos(par::AbstractArray{T, 1}, t::AbstractArray{T2, 1}) where {T<:Real, T2<:Real} = 
    doubleLog_template(par, t, doubleLog_Klos!)


export doubleLog_AG, doubleLog_Zhang, doubleLog_Beck, doubleLog_Elmore, doubleLog_Gu, doubleLog_Klos,
    doubleLog_AG!, doubleLog_Zhang!, doubleLog_Beck!, doubleLog_Elmore!, doubleLog_Gu!, doubleLog_Klos!
