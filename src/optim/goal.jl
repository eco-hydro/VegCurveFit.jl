function gof_RMSE(ypred::T, yobs::T, w::T=nothing) where {T<:AbstractArray{<:Real,1}}
  if w === nothing
    RMSE = sum((yobs - ypred) .^ 2.0)
  else
    RMSE = sum((yobs - ypred) .^ 2.0 .* w)
  end
  return (RMSE)
end

gof_RMSE(ypred::AbstractArray{<:Real,1}, input::input_struct) = gof_RMSE(ypred, input.y, input.w)


"""
    ypred = ones(input.y) * -0.1
    goal(par, doubleLog_Beck!, input)
    goal!(par, doubleLog_Beck!, input, ypred)

## Parameters
@param dot other parameters will be ignored
"""
function goal(par::T, FUN!::Function, input::input_struct) where {T<:Array{Float64,1}}
  ypred = ones(input.y) * 99.0
  goal!(par, FUN!, input, ypred)
end

# , ylu::Union{T,Nothing} = nothing
function goal!(par::T, FUN!::Function, input::input_struct, ypred::T) where {T<:Array{Float64,1}}
  FUN!(ypred, par, input.t)
  # check_ylu_w!(w, ypred, ylu)
  gof_RMSE(ypred, input)
end


# If ypred out of ylu, then set corresponding weights as `wmin`/
function check_ylu_w!(w::Union{T,Nothing}, ypred::T, ylu::Union{T,Nothing}) where {T<:Array{Float64,1}}
  if ylu !== nothing && w !== nothing
    n = length(ypred)
    @inbounds for i = 1:n
      if ypred[i] < ylu[1] || ypred[i] > ylu[2]
        w[i] = 0.0
      end
    end
  end
end


export goal!, goal
