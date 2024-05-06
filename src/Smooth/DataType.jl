# intermediate variables of whittaker smoother
import Parameters: @with_kw, @with_kw_noshow
export interm_whit


"""
    using struct object here will significantly improve SG performance 
    
```julia
par = interm_SG{Float32}(halfwin=5, d=3)
par = interm_SG(halfwin=5, d=3)
```
"""
@with_kw mutable struct interm_SG{FT}
  halfwin::Integer = 5
  frame::Integer = halfwin * 2 + 1
  d::Integer = 2

  S::AbstractArray{Int,2} = sgmat_S(halfwin, d)  # static array
  SMat::AbstractArray{FT,2} = Szeros(FT, frame, d + 1)   # [frame, d+1], temp variable, update everytime
  T::AbstractArray{FT,2} = Szeros(FT, d + 1, frame)   # [d+1  , frame]
  B::AbstractArray{FT,2} = Szeros(FT, frame, frame)  # [frame, frame]
  # n::Integer
  # y :: AbstractArray{FT, 1}
  # w :: AbstractArray{FT2, 1}
end


@with_kw mutable struct interm_whit{T}
  n::Int
  z::Vector{T} = zeros(T, n)

  # c: u1, d: v, e: u2
  d::Vector{T} = zeros(T, n) # v
  
  c::Vector{T} = zeros(T, n) # u1
  e::Vector{T} = zeros(T, n) # u2
  f::Vector{T} = zeros(T, n) # u2

  s0::Vector{T} = zeros(T, n)
  s1::Vector{T} = zeros(T, n)
  s2::Vector{T} = zeros(T, n)
end
