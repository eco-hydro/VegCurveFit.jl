function diffs!(xd::AbstractVector{T}, x::AbstractVector{T}, d=2) where {T<:Real}
  n = length(x)
  for k in 1:d
    if k == 1
      @inbounds for i in 1:n-k
        xd[i] = x[i+1] - x[i]
      end
    else
      @inbounds for i in 1:n-k
        xd[i] = xd[i+1] - xd[i]
      end
    end
  end
  @view(xd[1:n-d])
end

function diffs(x::AbstractVector{T}, d=2) where {T<:Real}
  d == 0 ? x : diffs(diff(x), d - 1)
end


function fidelity(y::AbstractVector, z::AbstractVector, w::AbstractVector)
  tol = 0.0f0
  @inbounds for i in eachindex(y)
    tol += w[i] * (y[i] - z[i])^2
  end
  log(tol)
end


function roughness!(zd::AbstractVector{T}, z::AbstractVector{T}, d=2) where {T<:Real}
  diffs!(zd, z, d)

  n = length(z)
  tol = T(0)
  @inbounds for i in 1:n-d
    tol += zd[i]^2
  end
  log(tol)
end

function roughness(z::AbstractVector{T}, d=2) where {T<:Real}
  zd = zeros(T, length(z))
  roughness!(zd, z, d)
end



"""
    lambda_vcurve(y::AbstractVector{T}, w::AbstractVector{T2};
        is_plot=false, lg_lambda_min=0.1, lg_lambda_max=3)

# Return
- `lambda`: optimal lambda
"""
function lambda_vcurve(y::AbstractVector{T}, w::AbstractVector{T2};
  interm=nothing,
  is_plot=false, lg_lambda_min=0.1, lg_lambda_max=3) where {
  T<:Real,T2<:Real}

  lg_lambdas = lg_lambda_min:0.1:lg_lambda_max
  n = length(lg_lambdas)
  fits = zeros(n)
  pens = zeros(n)

  # least of memory used
  # z = zeros(T, m)
  interm === nothing && (interm = interm_whit{Float32}(n=length(y)))
  zd = zeros(Float32, length(y))

  for i in 1:n
    lambda = 10^lg_lambdas[i]
    cve = whit2!(y, w, lambda, interm, include_cve=false)
    fits[i] = fidelity(y, interm.z, w)
    pens[i] = roughness!(zd, interm.z) # d=2
  end

  dfits = diff(fits)
  dpens = diff(pens)

  llastep = lg_lambdas[2] - lg_lambdas[1]
  v = @. sqrt(dfits^2 + dpens^2) / (log(10) * llastep)
  lamids = (lg_lambdas[2:end] + lg_lambdas[1:end-1]) / 2
  k = argmin(v)
  opt_lambda = 10^lamids[k]
  # z = whit2(y, lambda, w)
  if is_plot
    plot_lambda(y, w, lamids, v) |> display
  end
  opt_lambda
end


"""
    lambda_cv(y::AbstractVector{T}, w::AbstractVector{T2}; 
        interm = nothing, 
        is_plot=false, lg_lambda_min=0.1, lg_lambda_max=3)

# Return
- `lambda`: optimal lambda
"""
function lambda_cv(y::AbstractVector{T}, w::AbstractVector{T2};
  interm=nothing,
  is_plot=false, lg_lambda_min=0.1, lg_lambda_max=3) where {
  T<:Real,T2<:Real}

  lg_lambdas = lg_lambda_min:0.1:lg_lambda_max
  n = length(lg_lambdas)
  cvs = zeros(n)

  # least of memory used
  interm === nothing && (interm = interm_whit{Float32}(n=length(y)))

  # pens = zeros(n)
  for i in 1:n
    lambda = 10^lg_lambdas[i]
    # z, cvs[i] = whit2(y, w, lambda)
    cvs[i] = whit2!(y, w, lambda, interm; include_cve=true)
    # fits[i] = fidelity(y, z, w)
    # pens[i] = roughness(z, d)
  end

  # dfits = diff(fits)
  # dpens = diff(pens)
  # llastep = lg_lambdas[2] - lg_lambdas[1] 
  # v = @. sqrt(dfits^2 + dpens^2)/(log(10) * llastep)
  # lamids = (lg_lambdas[2:end] + lg_lambdas[1:end-1])/2
  k = argmin(cvs)
  opt_lambda = 10^lg_lambdas[k]
  if is_plot
    plot_lambda(y, w, lg_lambdas, cvs) |> display
  end
  opt_lambda
end


# x: lambdas candidates
function plot_lambda(y, w, lg_lambdas, cvs)
  k = argmin(cvs)
  opt_lambda = 10^lg_lambdas[k]

  p_v = plot(lg_lambdas, cvs, label="Generalized CV", frame=:box)
  scatter!(p_v, lg_lambdas, cvs, legend=false)
  scatter!(p_v, [lg_lambdas[k]], [cvs[k]],
    m=(10, :transparent, stroke(1, "red")),
    legend=false)
  vline!(p_v, [lg_lambdas[k]], color="red", linestyle=:dash)

  xlim = (0, length(y))
  p1 = plot(y, xlim=xlim, frame=:box)
  z_2, = whit2(y, w, 2.0)
  z_15, = whit2(y, w, 15.0)
  z_opt, = whit2(y, w, opt_lambda)
  plot!(p1, z_2, label="lambda = 2")
  plot!(p1, z_15, label="lambda = 15")
  plot!(p1, z_opt, label="lambda = $(round(opt_lambda, digits = 3))")
  plot(p_v, p1, layout=(1, 2), size=(700, 480))
end


export lambda_vcurve, lambda_cv, plot_lambda
