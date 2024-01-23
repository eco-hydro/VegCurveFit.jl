function doubleLog_Zhang!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  t0 = par[1]
  mn = par[2]
  mx = par[3]
  sos = par[4]
  rsp = par[5]
  eos = par[6]
  rau = par[7]

  if (t0 - sos <= 1 || t0 - eos >= -1)
    ypred .= 99.0
  end
  # pred <- (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
  #                     1 + exp( rau*(t[t >  t0] - eos))) + mn
  @inbounds for i = eachindex(t)
    if (t[i] <= t0)
      ypred[i] = mn + (mx - mn) / (1 + exp(-rsp * (t[i] - sos)))
    else
      ypred[i] = mn + (mx - mn) / (1 + exp(rau * (t[i] - eos)))
    end
  end
end

"""
doubleLog_AG
"""
function doubleLog_AG!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  t0 = par[1]
  mn = par[2]
  mx = par[3]
  rsp = par[4]
  a3 = par[5]
  rau = par[6]
  a5 = par[7]

  # pred <- mn + (mx - mn)*exp(- c( 
  #     ((t0 - t[t <= t0])*rsp) ^a3, 
  #     ((t[t >  t0] - t0)*rau) ^a5) )
  @inbounds for i = eachindex(t)
    if (t[i] <= t0)
      ypred[i] = mn + (mx - mn) * exp(-((t0 - t[i]) * rsp)^a3)
    else
      ypred[i] = mn + (mx - mn) * exp(-((t[i] - t0) * rau)^a5)
    end
  end
end

"""
doubleLog_Beck
"""
function doubleLog_Beck!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  # mn  = par[1]
  # mx  = par[2]
  sos = par[3]
  eos = par[5]
  rsp = par[4]
  rau = par[6]
  # pred <- mn + (mx - mn)*(
  #     1/(1 + exp(-rsp*(t - sos))) + 
  #     1/(1 + exp( rau*(t - eos))) - 1)
  @inbounds for i = eachindex(t)
    ypred[i] = par[1] + (par[2] - par[1]) * (
      1 / (1 + exp(-rsp * (t[i] - sos))) +
      1 / (1 + exp(rau * (t[i] - eos))) - 1)
  end
end

function doubleLog_Elmore!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  # mn  = par[1]
  # mx  = par[2]
  sos = par[3] # SOS
  rsp = par[4] # 1/rsp
  eos = par[5] # EOS
  rau = par[6] # 1/rau
  m7 = par[7]
  # pred <- mn + (mx - m7*t)*( 
  #     1/(1 + exp(-rsp*(t-sos))) - 
  #     1/(1 + exp(-rau*(t-eos))) )
  @inbounds for i = eachindex(t)
    ypred[i] = par[1] + (par[2] - m7 * t[i]) * (
      1 / (1 + exp(-rsp * (t[i] - sos))) -
      1 / (1 + exp(-rau * (t[i] - eos))))
  end
end

function doubleLog_Gu!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  y0 = par[1]
  a1 = par[2]
  a2 = par[3]
  sos = par[4]
  rsp = par[5]
  eos = par[6]
  rau = par[7]
  c1 = par[8]
  c2 = par[9]
  @inbounds for i = eachindex(t)
    # pred <- y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - 
    #              (a2/(1 + exp(-rau*(t - eos)))^c2)
    ypred[i] = y0 + (a1 / (1 + exp(-rsp * (t[i] - sos)))^c1) -
               (a2 / (1 + exp(-rau * (t[i] - eos)))^c2)
  end
end


function doubleLog_Klos!(ypred::AbstractVector{T}, par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real}
  a1 = par[1]
  a2 = par[2]
  b1 = par[3]
  b2 = par[4]
  c = par[5]
  B1 = par[6]
  B2 = par[7]
  m1 = par[8]
  m2 = par[9]
  q1 = par[10]
  q2 = par[11]
  v1 = par[12]
  v2 = par[13]
  @inbounds for i = eachindex(t)
    # pred <- (a1*t + b1) + (a2*t^2 + b2*t + c) * (
    #     1/(1 + q1 * exp(-B1 * (t - m1)))^v1 - 
    #     1/(1 + q2 * exp(-B2 * (t - m2)))^v2)
    ypred[i] = (a1 * t[i] + b1) + (a2 * t[i]^2 + b2 * t[i] + c) * (
      1 / (1 + q1 * exp(-B1 * (t[i] - m1)))^v1 -
      1 / (1 + q2 * exp(-B2 * (t[i] - m2)))^v2)
  end
end

## FUNCTIONS for USERS ---------------------------------------------------------

function doubleLog_template(par::AbstractVector{T}, t::AbstractArray{T2,1}, FUN!) where {T<:Real,T2<:Real}
  ypred = ones(T, length(t)) * 99.0
  FUN!(ypred, par, t)
  ypred
end

doubleLog_Zhang(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_Zhang!)

doubleLog_AG(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_AG!)

doubleLog_Beck(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_Beck!)

doubleLog_Elmore(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_Elmore!)

doubleLog_Gu(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_Gu!)
# c('y0', 'a1', 'a2', 'sos', 'rsp', 'eos', 'rau', 'c1', 'c2')

doubleLog_Klos(par::AbstractVector{T}, t::AbstractArray{T2,1}) where {T<:Real,T2<:Real} =
  doubleLog_template(par, t, doubleLog_Klos!)


export doubleLog_AG, doubleLog_Zhang, doubleLog_Beck, doubleLog_Elmore, doubleLog_Gu, doubleLog_Klos,
  doubleLog_AG!, doubleLog_Zhang!, doubleLog_Beck!, doubleLog_Elmore!, doubleLog_Gu!, doubleLog_Klos!
