using SparseArrays

speye(n) = SparseArrays.sparse(I, n, n)

function Base.diff(x::SparseMatrixCSC, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end


# TODO: bug here
function ddmat(x::AbstractVector, d::Integer=2)
  m = length(x)
  if d == 0
    return speye(m)
  else
    dx = x[(d+1):m] - x[1:(m-d)]
    V = spdiagm(0 => 1.0 ./ dx)
    return V * diff(ddmat(x, d - 1))
  end
end


function whit(y::AbstractVector, x::AbstractVector; lambda=2.0, d=2)
  m = length(y)
  E = speye(m)
  D = ddmat(x, d)

  M = E + lambda * D' * D
  C = cholesky(M, perm=1:m).L' # perm=1:3
  z = C \ (C' \ y)
  
  # H = inv(M)
  # h = diag(H)
  # r = (y - z) ./ (1 - h)
  # cve = sqrt(r' * r / m)

  # if m > 100
  #   n = 100
  #   E1 = speye(n)
  #   g = round.(Int, ((1:n) - 1) * (m - 1) / (n - 1) + 1)
    
  #   D1 = ddmat(x(g), d)

  #   lambda1 = lambda * (n / m)^(2 * d)
  #   H1 = inv(E1 + lambda1 * D1' * D1)
  #   h1 = diag(H1)
  #   u = zeros(m)
  #   k = floor(Int, m / 2)
  #   k1 = floor(Int, n / 2)

  #   u[k] = 1
  #   v = C \ (C' \ u)
  #   hk = v(k)
  #   f = round.(Int, ((1:m)' - 1) * (n - 1) / (m - 1) + 1)
  #   h = h1(f) * v(k) / h1(k1)
  # end
  z
end

export speye, ddmat, whit
