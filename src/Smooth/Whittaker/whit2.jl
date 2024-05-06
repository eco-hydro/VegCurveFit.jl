export whit2, whit2!, whit2_hat

"""
    whit2(y::AbstractVector{T}, w::AbstractVector{T2}, lambda::Float64; include_cve=true)

z, cve = whit2(y, w;lambda=2.0)
whit2(y, w; lambda)
whit2(y, w; lambda)
"""
function whit2(y::AbstractVector{T1}, w::AbstractVector{T2};
  lambda::Real, include_cve=true) where {T1<:Real,T2<:Real}

  interm = interm_whit{promote_type(T1, T2)}(; n=length(y))

  whit2!(y, w, lambda, interm; include_cve) # z, h, cve
end

"""
Second-order differences Whittaker-Henderson smoothing

LU decompose was used.

Ax = y
A = LU, LUx = y 
let `b = Ux`, `Lb = y`

# References
'Smoothing and interpolation with finite differences' [Eilers P. H. C, 1994]
(URL: http://dl.acm.org/citation.cfm?id=180916)
"""
function whit2!(y::AbstractVector{<:Real}, w::AbstractVector{<:Real}, lambda::Real, interm::interm_whit{FT};
  include_cve=true) where {FT<:Real}

  λ = FT(lambda)
  @unpack z, c, d, e = interm

  d[1] = w[1] + λ   # d是分母
  c[1] = -2λ / d[1] # 竖着数，第二项
  e[1] = λ / d[1]   # 竖着数，第三项
  z[1] = w[1] * y[1]

  d[2] = w[2] + 5λ - d[1] * c[1]^2 # 分母
  c[2] = (-4λ - d[1] * c[1] * e[1]) / d[2]
  e[2] = λ / d[2]
  z[2] = w[2] * y[2] - c[1] * z[1]

  # for (i = 2; i < m - 1; i++) 
  m = length(y)
  @inbounds @fastmath for i = 3:(m-2)
    i1 = i - 1
    i2 = i - 2
    d[i] = w[i] + 6λ - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
    c[i] = (-4λ - d[i1] * c[i1] * e[i1]) / d[i]
    e[i] = λ / d[i]
    z[i] = w[i] * y[i] - c[i1] * z[i1] - e[i2] * z[i2]
  end

  i = m - 1
  i1 = i - 1
  i2 = i - 2
  d[m-1] = w[m-1] + 5λ - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
  c[m-1] = (-2λ - d[i1] * c[i1] * e[i1]) / d[m-1]
  z[m-1] = w[m-1] * y[m-1] - c[i1] * z[i1] - e[i2] * z[i2]

  i = m
  i1 = i - 1
  i2 = i - 2
  d[m] = w[m] + λ - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
  z[m] = (w[m] * y[m] - c[i1] * z[i1] - e[i2] * z[i2]) / d[m]

  z[m-1] = z[m-1] / d[m-1] - c[m-1] * z[m]
  # for (i = m - 2; 0 <= i; i--)
  @inbounds @fastmath for i in (m-2):-1:1
    z[i] = z[i] / d[i] - c[i] * z[i+1] - e[i] * z[i+2]
  end

  cve::FT = ifelse(include_cve, whit2_hat(y, w, interm), FT(-999.0))
  z, cve
end


# according to hat and return the generalized cross validation
function whit2_hat(y::AbstractVector{<:Real}, w::AbstractVector{<:Real}, interm::interm_whit{FT}) where {FT<:Real}
  # # params: d, c, e, s0, s1, s2
  @unpack e, c, d = interm
  @unpack z, s0, s1, s2, n = interm

  # # Compute diagonal of inverse
  s0[n] = 1 / d[n]
  s0[n-1] = 1 / d[n-1] + c[n-1]^2 * s0[n]
  s1[n-1] = -c[n-1] * s0[n]
  
  @inbounds @fastmath for i = n-2:-1:1
    s1[i] = -c[i] * s0[i+1] - e[i] * s1[i+1]
    s2[i] = -c[i] * s1[i+1] - e[i] * s0[i+2]
    s0[i] = 1 / d[i] - c[i] * s1[i] - e[i] * s2[i]
  end

  tol = FT(0.0)
  wtol = FT(0.0)
  @inbounds for i = 1:n
    s0[i] *= w[i]
    r::FT = (y[i] - z[i]) / (1 - s0[i])
    tol += sum(r * r * w[i])
    wtol += w[i]
  end
  cve::FT = sqrt(tol / wtol) # sqrt(sum(r .* r) / n)
  cve
end
