using LinearAlgebra, SparseArrays
using Symbolics
import Symbolics: scalarize, variables

@variables λ

speye(n) = SparseArrays.sparse(I, n, n)

function Base.diff(x::SparseMatrixCSC, d::Integer=1)
  D = x[2:end, :] .- x[1:end-1, :]
  d >= 2 ? diff(D, d - 1) : D
end

function ddmat(x::AbstractVector, d::Integer=2)
  m = length(x)
  if d == 0
    return speye(m)
  else
    # dx = x[(d+1):m] - x[1:(m-d)] # bug may here
    return diff(ddmat(x, d - 1))
  end
end

function LU_decompose(A₁)
  n = size(A₁, 1)
  T = typeof(A₁)
  L = T(diagm(ones(n)))

  ## 徒手LU分解
  for i = 1:n-1
    r1 = A₁[i, :]
    # U[i, :] = r1
    for j = i+1:n
      f = A₁[j, i] / A₁[i, i]
      L[j, i] = f
      A₁[j, :] .= A₁[j, :] .- (f * r1)
      # 为啥要引入U，这样已经求解完成了
      # L[:, 1] = A₁[j, :]
    end
    # println("i = $i")
    # display(A₁)
  end
  (; L, U=A₁)
end

function diag_m(x)
  n = length(x)
  M = Matrix{Num}(undef, n, n)
  for i = 1:n, j = 1:n
    if i == j
      M[i, j] = x[i]
    else
      M[i, j] = 0.0
    end
  end
  M
end

# 代数余子式; algebraic complement
function complement(A::AbstractArray, i=1, j=1; verbose=false)
  i, j = j, i
  m, n = size(A)
  _A = A[setdiff(1:m, i), setdiff(1:n, j)]
  verbose && display(_A)
  (-1)^(i + j) * det(_A)
end

function complement(A::AbstractArray)
  R = similar(A)
  fill!(R, 0)
  m, n = size(A)
  for i = 1:m, j = 1:n
    R[i, j] = complement(A, i, j) # 代数余子式需要进行一次转置，才能得到A* = C'
  end
  R
end
