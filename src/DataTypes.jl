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
# TODO: 尚有提升空间
mutable struct input_struct{T<:Real}
  y::AbstractVector{T}
  # t::AbstractArray{<:Dates.Date, 1}
  t::AbstractArray{<:Union{Dates.Date,<:Real},1}
  w::AbstractVector{T}
end

function input_struct(y::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}
  input_struct{T}(y, t, ones(T, y))
end


mutable struct inputI_struct{T}
  y::AbstractVector{T}
  t::AbstractVector{<:Real}
  w::AbstractVector{T}
  ind::UnitRange{Int64}
  # tout::AbstractArray{<:Real, 1}
end


export input_struct, predictInput_struct, predictor_struct
