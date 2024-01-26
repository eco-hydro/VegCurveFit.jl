function wHANTS(y::AbstractVector, w, t; periodlen = 365, nf = 3, ignored...)
  n = length(y)
  # t = 1:n # TODO: 需要核对去除t的版本如何写

  ncol = min(2 * nf + 1, n)

  amp = zeros(nf + 1)
  phi = zeros(nf + 1)
  # z = zeros(n)

  ang = 2 * pi * (0:(periodlen-1)) / periodlen
  cs = cos(ang)
  sn = sin(ang)

  mat = zeros(n, ncol)
  mat[:, 1] = 1.0
  @inbounds for i = 1:nf
    I = 1 + mod(i * (t - 1), periodlen)
    mat[:, 2*i] = cs[I]
    mat[:, 2*i+1] = sn[I]
  end

  za = mat' * (w * y)
  A = (mat * w)' * mat # mat' * diag(w) * mat
  # % A = A + diag(ones(nr,1))*delta
  # % A(1,1) = A(1,1) - delta
  b = solve(A, za) # coefficients
  z = mat * b
  z = z[:, 1]

  amp[1] = b[1]
  phi[1] = 0.0
  i = 2:2:ncol
  ifr = (i + 2) / 2
  ra = b[i]
  rb = b[i+1]
  amp[ifr] = sqrt(ra * ra + rb * rb)

  dg = 180.0 / pi
  phase = atan2.(rb, ra) .* dg
  phase[phase<0] .= phase[phase<0] .+ 360
  phi[ifr] = phase
  
  z, amp, phi
end
