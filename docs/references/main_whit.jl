using LinearAlgebra, SparseArrays
using Symbolics

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
