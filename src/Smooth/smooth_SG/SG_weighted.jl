# using Whittaker2
using LinearAlgebra
using Parameters

as_SMatrix(mat::AbstractArray) = MMatrix{size(mat)...}(mat)
as_SMatrix(mat::SMatrix) = mat
as_SMatrix(mat::MMatrix) = mat

Szeros(FT::DataType, size...) = MMatrix{size...}(zeros(FT, size))
Szeros(size...) = Szeros(Float64, size...)


"""
    using struct object here will significantly improve SG performance 
    
```julia
par = SG_Param{Float32}(halfwin=5, d=3)
par = SG_Param(halfwin=5, d=3)
```
"""
@with_kw mutable struct SG_Param{FT}
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

function sgmat_wB!(w::AbstractArray{T1}, par::SG_Param) where {T1<:Real}
  # @show size(par.S) size(w)
  multiply_w_sqrt!(par.S, w, par.SMat)
  r = qr(par.SMat) # double
  par.T = as_SMatrix(r.R') \ par.S' # four times faster
  # par.T = r.R' \ par.S' # most time-consuing

  mul!(par.B, par.T', par.T)
  # par.B = par.T' * par.T;
  multiply_col!(par.B, w)
end

"""
    wSG(y::Array{T, 1}, w::Array{T2, 1}; halfwin=1, d=2)

weighted Savitzky Golay filter

# Arguments

- `check_wmin`: constrain the `w_min`. If not, it will lead to matrix division
  erorr.

# Examples

```julia
y = rand(100)
w = rand(100)
z1 = SG(y, halfwin = 5)
z2 = wSG(y, w, halfwin = 5)
```
"""
function wSG!(z::AbstractVector{FT}, y::AbstractVector{FT}, w::AbstractVector{FT};
  halfwin=1, d=2, par=nothing, check_wmin=false) where {FT<:Real}

  

  n = length(y)
  if par === nothing
    par = SG_Param{FT}(halfwin=halfwin, d=d)
  end
  d = par.d
  frame = par.frame
  halfwin = par.halfwin

  if check_wmin
    w = check_wmin(w)
  end

  if (sum(w) == n)
    z .= SG(y; halfwin=par.halfwin, d)
    return z
  end

  # y_head = @views(par.B[1:halfwin+1, :] * y[1:frame])[:, 1]
  sgmat_wB!(w[1:frame], par)
  @inbounds @views for i = 1:halfwin+1
    z[i] = dot(par.B[i, :], y[1:frame])
  end

  # y_mid = zeros(eltype(y), n - frame - 1) # vector
  @inbounds @views for i = 1:n-frame-1
    sgmat_wB!(w[i+1:i+frame], par)
    k = i + halfwin + 1
    z[k] = dot(par.B[halfwin+1, :], y[i+1:i+frame])
  end

  # y_tail = @views(par.B[halfwin+1:frame, :] * y[n-frame+1:n])[:, 1]
  sgmat_wB!(w[n-frame+1:n], par)
  @inbounds @views for i = halfwin+1:frame
    k = i - frame + n
    z[k] = dot(par.B[i, :], y[n-frame+1:n])
  end
end

function wSG(y::AbstractVector{FT}, w::AbstractVector{FT}; kw...) where {FT<:Real}
  z = zeros(FT, size(y))
  wSG!(z, y, w; kw...)
  z
end

wSG(y::AbstractVector{FT}; kw...) where {FT<:Real} = wSG(y, ones(FT, size(y)); kw...)


# weighted Savitzky Golay filter
function wSG_low(y::AbstractVector{T}, w::AbstractArray{T2,1}; halfwin=1, d=2) where {
  T<:Real,T2<:Real}
  # constrain the w_min, unless it will lead to matrix division erorr
  # w = deepcopy(w)
  # wmin = 1e-4
  # w[w.<wmin] .= wmin

  n = length(y)
  frame = halfwin * 2 + 1
  if (sum(w) == n)
    return SG(y; halfwin=halfwin, d=d)
  end

  S = sgmat_S(halfwin, d)
  smat = ones(T2, size(S))

  B = sgmat_wB(S, w[1:frame], smat)
  y_head = @views(B[1:halfwin+1, :] * y[1:frame])[:, 1]

  y_mid = zeros(T, n - frame - 1)
  @inbounds for i = 1:n-frame-1
    B = sgmat_wB(S, w[i+1:i+frame], smat)
    y_mid[i] = dot(@view(B[halfwin+1, :]), @view y[i+1:i+frame])
  end

  B = sgmat_wB(S, w[n-frame+1:n], smat)
  y_tail = @views(B[halfwin+1:frame, :] * y[n-frame+1:n])[:, 1]
  [y_head; y_mid; y_tail]
end


export Szeros, as_SMatrix,
  SG_Param, wSG, wSG_low, wSG_old;
