
mutable struct param_struct
  mx
  mn
  ampl
  sos
  eos
  doy_peak
  deltaT
  deltaY
  half
  t1
  t2
  k
end

# method, par, tout, ylu
mutable struct predictInput_struct
  method::Vector{String} # vector of string
  par::Vector{Vector{Float64}}   # vector of model parameters
  tout::UnitRange{Int64}
  ylu::Vector{Float64}
end

mutable struct predictor_struct
  # method::Vector{String} # vector of string
  opt::Vector{Dict{String,Any}}   # vector of model parameters
  tout::UnitRange{Int64}
  ylu::Vector{Float64}
end



# input of check_input object
mutable struct input_struct
  y::AbstractArray{<:Real,1}
  # t::AbstractArray{<:Dates.Date, 1}
  t::AbstractArray{<:Union{Dates.Date,<:Real},1}
  w::AbstractArray{<:Real,1}
end

function input_struct(y::AbstractArray{T,1}, t::AbstractArray{T,1}) where {T<:Real}
  input_struct(y, t, ones(y))
end

mutable struct inputI_struct
  y::AbstractArray{<:Real,1}
  t::AbstractArray{<:Real,1}
  w::AbstractArray{<:Real,1}
  ind::UnitRange{Int64}
  # tout::AbstractArray{<:Real, 1}
end


export input_struct, predictInput_struct, predictor_struct
