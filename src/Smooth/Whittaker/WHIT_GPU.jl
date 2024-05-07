using CUDA

function WHIT_GPU(y::AbstractVector, w::AbstractVector, x::AbstractVector;
  lambda=2.0, d=2, ignore...)
  # Move data to GPU
  y = CUDA.CuArray(y)
  w = CUDA.CuArray(w)
  x = CUDA.CuArray(x)

  n = length(y)
  D = ddmat(x, d)

  W = SparseArrays.spdiagm(w)
  A = W + lambda * D' * D

  # L = cholesky(A).L # Matrix
  L = cholesky(A, perm=1:n).L # sparse
  z = L' \ (L \ (w .* y))
  z
end
# export WHIT, speye, diff, ddmat
