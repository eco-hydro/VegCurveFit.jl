"""
Type-2 kurtosis and skewness
"""
function kurtosis(x::AbstractVector{T}) where {T<:Real}
  n = length(x)
  x = x .- mean(x)
  r = n * sum(x .^ 4) / (sum(x .^ 2)^2)

  ((n + 1) * (r - 3) + 6) * (n - 1) / ((n - 2) * (n - 3)) # type_2
  # r * (1 - 1/n)^2 - 3 # type_3
end

function skewness(x::AbstractVector{T}) where {T<:Real}
  n = length(x)
  x = x .- mean(x)
  r = sqrt(n) * sum(x .^ 3) / (sum(x .^ 2)^(3 / 2))
  r * sqrt(n * (n - 1)) / (n - 2) # type2
  # r * ((1 - 1/n))^(3/2) # type3
end

"""
Init the Whittaker parameter lambda 

- coef: coefficients in the order of interception, mean, sd, skewness and kurtosis, 
    e.g. [0.9809, 0.7247, -2.6752, -0.3854, -0.0604].

# Examples
```jldoctest
y = rand(100)
coef = [0.9809, 0.7247, -2.6752, -0.3854, -0.0604]; # Terra EVI
lambda_init(y, coef)
```
"""
function lambda_init(x::AbstractVector{T}, 
  coef::AbstractVector{Float64}=[0.9809, 0.7247, -2.6752, -0.3854, -0.0604]) where {T<:Real}

  x_mean = mean(x)
  x_sd = std(x)

  # cv = x_sd / x_mean
  skew = skewness(x)
  kur = kurtosis(x)

  lambda = coef[1] + coef[2] * x_mean + coef[3] * x_sd + coef[4] * skew + coef[5] * kur 
  # lambda = 0.9809 + 0.7247 * x_mean - 2.6752 * x_sd - 0.3854 * skew - 0.0604 * kur
  10^lambda
end


export kurtosis, skewness, lambda_init

# x = rand(100)
# x[3] = missing
# lambda_init(x)
# R"phenofit::lambda_init($x)"

# std(x)
# R"sd($x)"
# R"phenofit::kurtosis($x)"
# R"phenofit::skewness($x)"
