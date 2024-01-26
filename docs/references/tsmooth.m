function [x, score] = tsmooth(y, sig, J)

n = length(y); nc = ceil(n/2);
elim = 2 *(1 - sig); flim =(1 - sig)/(1 + sig); lam = 4 * sig^4/(1 - sig^2);
N = ceil(1 - J/ log 10(flim)); glim =(1 - sig^2)/(4 * sig^3 *(2 - sig^2));
e = zeros(1, N + 1); f = zeros(1, N + 2); x = zeros(1, n);
a1 = 1 + lam; a2 = 5 + lam; a3 = 6 + lam;

if N > nc
  error('sig too small, use smooth instead')
end

%Factor the coefficient matrix and solve the first triangular system
d = a1; f(1) = 1/d; x(1) = f(1) * lam * y(1); mu = 2; e(1) = mu * f(1);
d = a2 - mu * e(1); f(2) = 1/d; x(2) = f(2) *(lam * y(2) + mu * x(1)); mu = 4 - e(1); e(2) = mu * f(2);
for j = 3 : N
  m1 = j - 1; m2 = j - 2;
  d = a3 - mu * e(m1) - f(m2);
  f(j) = 1/d;
  x(j) = f(j) *(lam * y(j) + mu * x(m1) - x(m2));
  mu = 4 - e(m1);
  e(j) = mu * f(j);
end
mu = 4 - elim;
for j = N + 1 : n - 2
  x(j) = flim * (lam * y(j) + mu * x(j - 1) - x(j - 2));
end

d = a2 - mu * elim - flim; f(N + 1) = 1/d;
x(n - 1) = f(N + 1) *(lam * y(n - 1) + mu * x(n - 2) - x(n - 3));
mu = 2 - elim; e(N + 1) = mu * f(N + 1);
d = a1 - mu * e(N + 1) - flim; f(N + 2) = 1/d;
x(n) = f(N + 2) *(lam * y(n) + mu * x(n - 1) - x(n - 2));
%Solve the second triangular system and find avg squared error
sq =(y(n) - x(n))^2;
x(n - 1) = x(n - 1) + e(N + 1) * x(n);
sq = sq + (y(n - 1) - x(n - 1))^2;
for j = n - 2 : -1 : N + 1
  x(j) = x(j) + elim * x(j + 1) - flim * x(j + 2);
  sq = sq + (y(j) - x(j))^2;
end
for j = N : -1 : 1
  x(j) = x(j) + e(j) * x(j + 1) - f(j) * x(j + 2);
  sq = sq + (y(j) - x(j))^2;
end
sq = sq/n;

%Compute GCV score
g2 = f(N + 2); tr = g2; h = e(N + 1) * g2;
g1 = f(N + 1) + e(N + 1) * h; tr = tr + g1;
for k = n - 2 : -1 : n - N + 1
  q = elim * h - flim * g2;
  h = elim * g1 - flim * h; g2 = g1;
  g1 = flim + elim * h - flim * q;
  tr = tr + g1;
end

tr = tr + (nc - N) * glim;
tr =(2 * tr - rem(n, 2) * glim) * lam/n;
score = sq/(1 - tr)^2;

end
