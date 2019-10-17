"""
doubleLog_Beck
"""
function doubleLog_Beck(par::T, t::T)::T where T <: Array{Float64,1}
    mn  = par[1]
    mx  = par[2]
    sos = par[3]
    rsp = par[4]
    eos = par[5]
    rau = par[6]
    # if (sos >= eos) return(rep(9999, length(t)))
    if (eos < sos); return(ones(length(t)) * 99.0); end

    pred = mn .+ (mx - mn) * (1 ./ (1 .+ exp.(-rsp * (t .- sos))) + 
        1 ./ (1 .+ exp.(rau * (t .- eos))) .- 1)
end

function doubleLog_Beck!(ypred::T, par::T, t::T) where T <: Array{Float64,1}
    sos = par[3]
    eos = par[5]
    
    n = length(t)
    # if (sos >= eos) return(rep(9999, length(t)))
    if (eos < sos)
        @inbounds for i = 1:n
            ypred[i] = 99.0
        end
        return
    else
        # mn  = par[1]
        # mx  = par[2]
        rsp = par[4]
        rau = par[6]
        
        @inbounds for i = 1:n
            ypred[i] = par[1] + (par[2] - par[1]) * (1 / (1 + exp(-rsp * (t[i] - sos))) + 
                1 / (1 + exp(rau * (t[i] - eos))) - 1)
        end
    end
end
