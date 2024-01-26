# using Whittaker2
using LinearAlgebra
using Parameters

as_SMatrix(mat::AbstractArray) = MMatrix{size(mat)...}(mat)
as_SMatrix(mat::SMatrix) = mat
as_SMatrix(mat::MMatrix) = mat

Szeros(FT::DataType, size...) = MMatrix{size...}(zeros(FT, size))
Szeros(size...) = Szeros(Float64, size...)


function sgmat_wB!(w::AbstractArray{T1}, interm::interm_SG) where {T1<:Real}
  # @show size(interm.S) size(w)
  multiply_w_sqrt!(interm.S, w, interm.SMat)
  r = qr(interm.SMat) # double
  interm.T = as_SMatrix(r.R') \ interm.S' # four times faster
  # interm.T = r.R' \ interm.S' # most time-consuing

  mul!(interm.B, interm.T', interm.T)
  # interm.B = interm.T' * interm.T;
  multiply_col!(interm.B, w)
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
  halfwin=1, d=2, check_wmin=false, 
  interm::interm_SG=interm_SG{FT}(; halfwin, d), ignored...) where {FT<:Real}

  n = length(y)
  @unpack d, frame, halfwin = interm

  check_wmin && (w = _check_wmin(w))

  if (sum(w) == n)
    z .= SG(y; halfwin=interm.halfwin, d)
    return z
  end

  # y_head = @views(interm.B[1:halfwin+1, :] * y[1:frame])[:, 1]
  sgmat_wB!(w[1:frame], interm)
  @inbounds @views for i = 1:halfwin+1
    z[i] = dot(interm.B[i, :], y[1:frame])
  end

  # y_mid = zeros(eltype(y), n - frame - 1) # vector
  @inbounds @views for i = 1:n-frame-1
    sgmat_wB!(w[i+1:i+frame], interm)
    k = i + halfwin + 1
    z[k] = dot(interm.B[halfwin+1, :], y[i+1:i+frame])
  end

  # y_tail = @views(interm.B[halfwin+1:frame, :] * y[n-frame+1:n])[:, 1]
  sgmat_wB!(w[n-frame+1:n], interm)
  @inbounds @views for i = halfwin+1:frame
    k = i - frame + n
    z[k] = dot(interm.B[i, :], y[n-frame+1:n])
  end
  z
end


wSG_low(y::AbstractVector{FT}; kw...) where {FT<:Real} =
  wSG_low(y, ones(FT, size(y)); kw...)


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
  interm_SG, wSG, wSG_low, wSG_old;
