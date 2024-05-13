using SparseArrays

function Base.diff(x::AbstractMatrix, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diffn(D, d - 1) : D
end

function Base.diff(x::SparseMatrixCSC, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end

# ddmat(x::AbstractVector, d::Integer=2) = diff(spdiagm(x), d)

# 相比于原版进行了微调
function ddmat(x::AbstractVector, d::Integer=2)
  n = length(x)
  if d == 0
    return sparse(I, n, n)
  else
    dx = ( x[d+1:n] - x[1:(n-d)] ) ./ d
    V = spdiagm( 1 ./ dx)
    return V * diff(ddmat(x, d - 1))
  end
end

function WHIT(y::AbstractVector, w::AbstractVector; kw...)
  WHIT(y, w, 1:length(y); kw...)
end

function WHIT(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  lambda=2.0, d=2, include_cve=true)
  n = length(y)
  D = ddmat(x, d)

  W = SparseArrays.spdiagm(w)
  A = W + lambda * D' * D

  # L = cholesky(A).L # Matrix
  L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w.*y))
  # display(sparse(L)[1:7, 1:7])

  cve = -999.0
  if include_cve
    inv_L = Matrix(sparse(L))^-1 # 这一步可能消耗了较多的时间
    H = inv_L' * inv_L * W # 借力cholesky

    h = diag(H)
    r = @. (y - z) / (1 - h)
    cve = sqrt(sum(r .* r .* w) / sum(w))
  end
  z, cve
end

export WHIT, speye, diff, ddmat

