# 一种`带状矩阵对角阵`的快速算法
# 用于`Whittaker smoother`求解
# 
# Dongdong Kong, CUG, 2024-05-07

function def_U(c, e, f)
  n = length(c)
  U = Matrix{Num}(undef, n, n)
  fill!(U, 0)
  for i = 1:n
    U[i, i] = 1
    i + 1 <= n && (U[i, i+1] = c[i])
    i + 2 <= n && (U[i, i+2] = e[i])
    i + 3 <= n && (U[i, i+3] = f[i])
  end
  U
end

function def_U(c, e)
  n = length(c)
  U = Matrix{Num}(undef, n, n)
  fill!(U, 0)
  for i = 1:n
    U[i, i] = 1
    i + 1 <= n && (U[i, i+1] = c[i])
    i + 2 <= n && (U[i, i+2] = e[i])
  end
  U
end

"""
- `U`: [
  1 c₁ e₁ f₁ 0
  0 1  c₂ e₂ f₂
  0 0  1  c₃ e₃
  0 0  0  1  c₄
  0 0  0  0  1
]
"""
function cal_diag(U2, d; m=3)
  # S2: [n, m+1]
  # U2: [n, m]
  n = length(d)
  S = variables(:S, 1:n, 1:m+1) # m=2,3个临时变量已足够
  fill!(S, 0)
  S[n, 1] = 1 / d[n]

  for i = n-1:-1:1
    S[i, 1] = 1 / d[i]
    for l = 1:min(m, n - i)
      S[i, 1+l] = 0
      for k = 1:min(n - i, m)
        # if k <= l
        #   S[i, 1+l] -= U[i, i+k] * S[i+k, l-k+1]
        # else
        #   S[i, 1+l] -= U[i, i+k] * S[i+l, k-l+1]
        # end
        _i, _j = k <= l ? (i + k, l - k + 1) : (i + l, k - l + 1)
        S[i, 1+l] -= U2[i, k] * S[_i, _j]
      end
      S[i, 1] -= U2[i, l] * S[i, 1+l]
    end
  end
  S[:, 1]
end

function cal_diag_full(U, d; m=3)
  # B = (U' * D * U)                        # Hutchinson 1985, Eq. 3.1
  # B^(-1) = B * U^(-1)' + (1 - U) * B^(-1) # Hutchinson 1985, Eq. 3.3
  n = length(d)
  b = variables(:b, 1:n, 1:n) # elements of B^-1, 对称矩阵，因此只求一半的元素即可
  b[n, n] = 1 / d[n]
  # b[n-1, n] = -U[n-1,n]*b[n,n]
  # b[n-1, n-1] = d[n-1] - U[n-1,n]*b[n-1,n]  

  for i = n-1:-1:1
    b[i, i] = 1 / d[i]

    for l = 1:min(m, n - i)
      b[i, i+l] = 0
      for k = 1:min(n - i, m)
        _i, _j = k <= l ? (i + k, i + l) : (i + l, i + k)
        b[i, i+l] -= U[i, i+k] * b[_i, _j]
        # if k <= l
        #   b[i, i+l] -= U[i, i+k] * b[i+k, i+l]
        # else
        #   b[i, i+l] -= U[i, i+k] * b[i+l, i+k]
        # end
      end
      b[i, i] -= U[i, i+l] * b[i, i+l]
    end
  end
  diag(b)
end
