
function HANTS(y, t, w; periodlen = 365, nf = 3)
  n = length(y)
  ncol = min(2 * nf + 1, n)

  amp = zeros(nf + 1)
  phi = zeros(nf + 1)
  # yz = zeros(n)

  ang = 2 * pi * (0:(periodlen-1)) / periodlen
  cs = cos(ang)
  sn = sin(ang)

  mat = zeros(n, ncol)
  mat[:, 1] = 1.0
  for i = 1:nf
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


function smooth_HANTS(y, qc, date;
  FUN = HANTS,
  wFUN=wBisquare,
  
  nptperyear=46,
  iters=3,
  param = (nf=3, periodlen=365),
  param_wFUN=(;),
  wmin=0.2, wmid=0.5, wmax=1.0,
  alpha=0.02,
  ignored...)
  
  w, QC_flag = qc_FparLai(qc; wmin=wmin, wmid=wmid, wmax=wmax)
  ylu, wc = get_ylu(y, w; wmin=wmin, wmid=wmid, wmax=wmax, alpha=alpha)

  data = DataFrame(; date, y, w, QC_flag)

  # w1 = ones(size(y));
  res = map(i -> begin
      yfit = FUN(y, w; param...)
      clamp!(yfit, ylu[1], Inf) # constrain in the range of ylu

      w = wFUN(y, yfit, w; iter=i, nptperyear=nptperyear, wmin=0.05, param_wFUN...)
      dfit = DataFrame(; date, z=yfit, w)
    end, 1:iters)

  predict = melt_list(res, iter=1:iters)
  Dict("data" => data, "predict" => predict,
    "param" => param)
end
