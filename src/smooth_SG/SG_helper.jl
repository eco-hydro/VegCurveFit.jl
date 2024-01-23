using StaticArrays


function sgmat_S(halfwin::Int=1, d::Int=2)
  frame = 2 * halfwin + 1
  mat = zeros(Int, frame, d + 1)
  # mat = zeros(SMatrix{frame, d+1})    
  for i = 0:frame-1, j = 0:d
    mat[i+1, j+1] = (i - halfwin)^j # fix solaris error
  end
  mat
end

# Update Smat matrix in high efficiency way, reuse smat
# multiply_row
"""
  multiply_w_sqrt(S::AbstractMatrix, w::AbstractVector)

```julia
S = rand(11, 3)
w = rand(11)

multiply_w_sqrt(S, w)

## second version
S = MMatrix{11, 3}(rand(11, 3))
w = rand(11)
multiply_w_sqrt(S, w)
```
"""
function multiply_w_sqrt(S::AbstractMatrix, w::AbstractVector)
  # repeat(sqrt.(w), 1, size(S, 2)) .* S
  Snew = zeros(eltype(w), size(S))
  multiply_w_sqrt!(S, w, Snew)
  Snew
end

function multiply_w_sqrt(S::MMatrix, w::AbstractVector)
  _size = size(S)
  Snew = MMatrix{_size...}(zeros(eltype(S), _size))
  multiply_w_sqrt!(S, w, Snew)
  Snew
end

function multiply_w_sqrt!(S::AbstractMatrix{T}, w::AbstractVector{FT}, Snew::AbstractMatrix{FT}) where {FT<:Real, T<:Real}
  nrow, ncol = size(S)
  @inbounds for i in 1:nrow, j in 1:ncol
    Snew[i, j] = S[i, j] * sqrt(w[i])
  end
end

function check_wmin(w::AbstractVecOrMat; wmin=1e-4)
  w = deepcopy(w)
  w[w.<wmin] .= wmin
  w
end


export multiply_w_sqrt, check_wmin;
