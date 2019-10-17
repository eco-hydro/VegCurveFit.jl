function gof_RMSE(yobs::T, ypred::T, w::Union{T,Nothing} = nothing)::Float64 where T<:Array{Float64,1}

    if w === nothing
        RMSE = sum((yobs - ypred).^ 2.0)
    else
        RMSE = sum((yobs - ypred).^ 2.0 .* w)
    end
    return(RMSE)
end


"""
@param dot other parameters will be ignored
"""
function goal(par::T, FUN::Function, y::T, t::T, 
    w::Union{T,Nothing} = nothing, ylu::Union{T,Nothing} = nothing) where T<:Array{Float64,1}

    ypred = FUN(par, t)
    check_ylu_w!(w, ypred, ylu)
    gof_RMSE(y, ypred, w)
end

function goal!(par::T, FUN::Function, y::T, t::T, ypred::T, 
    w::Union{T,Nothing} = nothing, ylu::Union{T,Nothing} = nothing) where T<:Array{Float64,1}

    FUN(ypred, par, t)    
    check_ylu_w!(w, ypred, ylu)
    gof_RMSE(y, ypred, w)
end

# If ypred out of ylu, then set corresponding weights as `wmin`/
function check_ylu_w!(w::Union{T, Nothing}, ypred::T, ylu::Union{T, Nothing}) where T<:Array{Float64,1}

    if ylu !== nothing && w !== nothing
        n = length(ypred)
        @inbounds for i = 1:n
            if ypred[i] < ylu[1] || ypred[i] > ylu[2]
                w[i] = 0.0
            end
        end
    end
end
