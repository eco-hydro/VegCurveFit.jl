multiply_col!(S::AbstractArray{T1,2}, w::AbstractArray{T}) where {
  T1<:Real,T<:Real} = begin
  nrow, ncol = size(S)
  @inbounds for i in 1:nrow, j in 1:ncol
    S[i, j] = S[i, j] * w[j] # note is w[j], not 
  end
end

"""
    B matrix of Savitzky Golay

sgmat_B(S::Array{Int, 2})   
sgmat_wB(S::AbstractArray{Int, 2}, w)
"""
function sgmat_B(S::Array{Int,2})
  r = qr(S)
  T = r.R' \ S'
  B = T' * T
  return B
end

# B matrix of weighted Savitzky Golay
function sgmat_wB(S::AbstractArray{Int,2}, w)
  r = qr(multiply_w_sqrt(S, w))
  # T = as_SMatrix(r.R') \ S'; # four times faster
  T = r.R' \ S' # most time-consuing
  B = T' * T
  B = repeat(w, 1, size(B, 2))' .* B
  return B
end

function sgmat_wB(S::AbstractArray{Int,2}, w::AbstractArray{T1}, smat::AbstractArray{T1,2}) where {T1<:Real}
  multiply_w_sqrt!(S, w, smat)
  r = qr(smat) # double
  # T = as_SMatrix(r.R') \ S'; # four times faster
  T = r.R' \ S' # most time-consuing

  B = T' * T
  multiply_col!(B, w)
  return B
end

"""
    SG(y::Array{T, 1}; halfwin=1, d=2)   
    SG(y::Array{T, 1}, w::Array{T2, 1}; halfwin=1, d=2)

weighted Savitzky Golay filter

# Examples
y = rand(100)
w = rand(100)
z1 = SG(y, halfwin = 5)
z2 = wSG(y, w, halfwin = 5)
"""
# Savitzky Golay filter
function SG(y::Array{T,1}; halfwin=1, d=2) where {T<:Real}
  frame = halfwin * 2 + 1

  S = sgmat_S(halfwin, d)
  B = sgmat_B(S)
  y_head = @view(B[1:halfwin+1, :]) * @view y[1:frame]

  n = length(y)
  y_mid = zeros(T, n - frame - 1, 1)
  @inbounds for i = 1:n-frame-1
    y_mid[i] = dot(@view(B[halfwin+1, :]), @view y[i+1:i+frame])
    # y_mid[i] = dot(B[halfwin+1, :], y[i:i+frame-1]);
  end
  y_tail = @view(B[halfwin+1:frame, :]) * @view y[n-frame+1:n]
  [y_head; y_mid; y_tail][:, 1]
end

export sgmat_S, sgmat_B, sgmat_wB, SG;
