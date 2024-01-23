# 定义一个中间变量盛放中间变量
# intermediate variables of whittaker smoother
import Parameters: @with_kw, @with_kw_noshow

@with_kw mutable struct interm_whit{T}
  n::Int
  z::Vector{T} = zeros(T, n)

  # c: u1, d: v, e: u2
  c::Vector{T} = zeros(T, n) # u1
  d::Vector{T} = zeros(T, n) # v
  e::Vector{T} = zeros(T, n) # u2

  s0::Vector{T} = zeros(T, n)
  s1::Vector{T} = zeros(T, n)
  s2::Vector{T} = zeros(T, n)
end


export interm_whit



"""
    whit2(y::AbstractVector{T}, w::AbstractVector{T2}, lambda::Float64; include_cve=true)

Second-order differences Whittaker-Henderson smoothing

z, cve = whit2(y, w, lambda)
whit2(y, w, lambda, z)
whit2(y, w, lambda, z, c, d, e)

# Citation
'Smoothing and interpolation with finite differences' [Eilers P. H. C, 1994]
(URL: http://dl.acm.org/citation.cfm?id=180916)
"""
function whit2(y::AbstractVector{<:Real}, w::AbstractVector{FT}, lamb::Real; include_cve=true) where {
  FT<:Real}

  # FT = Float32
  interm = interm_whit{FT}(; n=length(y))

  cve = whit2!(y, w, lamb, interm; include_cve)
  interm.z, cve
end

# function whit2!(y::AbstractVector{T}, w::AbstractVector{T2}, lambda::Float64, z::AbstractVector{T}; include_cve=true) where {
#   T<:Real,T2<:Real}
#   n = length(y)
#   interm = interm_whit{T}(; n)
#   cve = whit2!(y, w, lambda, z, interm; include_cve)
#   cve
# end

function whit2!(y::AbstractVector{<:Real}, w::AbstractVector{FT}, lamb::Real, interm::interm_whit{FT};
  include_cve=true) where {FT<:Real}

  lambda::FT = FT(lamb)
  @unpack z, c, d, e = interm
  # z = interm.z
  # c = interm.c
  # d = interm.d
  # e = interm.e
  d[1] = w[1] + lambda
  c[1] = -2 * lambda / d[1]
  e[1] = lambda / d[1]
  z[1] = w[1] * y[1]
  d[2] = w[2] + 5 * lambda - d[1] * c[1] * c[1]
  c[2] = (-4 * lambda - d[1] * c[1] * e[1]) / d[2]
  e[2] = lambda / d[2]
  z[2] = w[2] * y[2] - c[1] * z[1]

  # for (i = 2; i < m - 1; i++) 
  m = length(y)
  @inbounds @fastmath for i = 3:(m-1)
    i1 = i - 1
    i2 = i - 2
    d[i] = w[i] + 6 * lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
    c[i] = (-4 * lambda - d[i1] * c[i1] * e[i1]) / d[i]
    e[i] = lambda / d[i]
    z[i] = w[i] * y[i] - c[i1] * z[i1] - e[i2] * z[i2]
  end

  i = m - 1
  i1 = i - 1
  i2 = i - 2
  d[m-1] = w[m-1] + 5 * lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
  c[m-1] = (-2 * lambda - d[i1] * c[i1] * e[i1]) / d[m-1]
  z[m-1] = w[m-1] * y[m-1] - c[i1] * z[i1] - e[i2] * z[i2]

  i = m
  i1 = i - 1
  i2 = i - 2
  d[m] = w[m] + lambda - c[i1] * c[i1] * d[i1] - e[i2] * e[i2] * d[i2]
  z[m] = (w[m] * y[m] - c[i1] * z[i1] - e[i2] * z[i2]) / d[m]
  z[m-1] = z[m-1] / d[m-1] - c[m-1] * z[m]

  # for (i = m - 2; 0 <= i; i--)
  @inbounds @fastmath for i in (m-2):-1:1
    z[i] = z[i] / d[i] - c[i] * z[i+1] - e[i] * z[i+2]
  end

  # cve = -990.0
  cve::FT = ifelse(include_cve, whit2_hat(y, w, interm), FT(-999.0))
  cve
end


# according to hat and return the generalized cross validation
function whit2_hat(y::AbstractVector{<:Real}, w::AbstractVector{FT}, interm::interm_whit{FT}) where {FT<:Real}
  # c: u1, d: v, e: u2
  # # params: v, u1, u2, s0, s1, s2
  v = interm.d
  u1 = interm.c
  u2 = interm.e
  @unpack z, s0, s1, s2, n = interm
  # z = interm.z
  # s0 = interm.s0
  # s1 = interm.s1
  # s2 = interm.s2  

  # # Compute diagonal of inverse
  @inbounds @fastmath for i = n:-1:1
    i1 = i + 1
    i2 = i + 2
    s0[i] = 1 / v[i]
    if (i < n)
      s1[i] = -u1[i] * s0[i1]
      s0[i] = 1 / v[i] - u1[i] * s1[i]
    end
    if (i < n - 1)
      s1[i] = -u1[i] * s0[i1] - u2[i] * s1[i1]
      s2[i] = -u1[i] * s1[i1] - u2[i] * s0[i2]
      s0[i] = 1 / v[i] - u1[i] * s1[i] - u2[i] * s2[i]
    end
  end

  tol = FT(0.0)
  @inbounds for i = 1:n
    # r = @. (y - z) * w / (1 - s0) # 这一步生成了一个新的变量
    r::FT = (y[i] - z[i]) * w[i] / (1 - s0[i])
    tol += sum(r * r)
  end

  cve::FT = sqrt(tol / n) # sqrt(sum(r .* r) / n)
  cve
end

export whit2, whit2!,
  whit2_cpp,
  whit2_hat
