function [x, score] = smooth(y, lam)

n = length(y); nc = ceil(n/2);
e = zeros(1, n - 1); f = zeros(1, n); x = zeros(1, n);
% lam = 4 * sig^4/(1 - sig^2);

a1 = 1 + lam; 
a2 = 5 + lam; 
a3 = 6 + lam;

%Factor the coefficient matrix and solve the first triangular system
d = a1; 
f(1) = 1/d;
x(1) = f(1) * lam * y(1);

mu = 2; 
e(1) = mu * f(1);
d = a2 - mu * e(1);
f(2) = 1/d;
x(2) = f(2) *(lam * y(2) + mu * x(1));

mu = 4 - e(1);
e(2) = mu * f(2);

for j = 3 : n - 2
  m1 = j - 1;
  m2 = j - 2;
  d = a3 - mu * e(m1) - f(m2);
  f(j) = 1/d;
  x(j) = f(j) *(lam * y(j) + mu * x(m1) - x(m2));
  mu = 4 - e(m1);
  e(j) = mu * f(j);
end

d = a2 - mu * e(n - 2) - f(n - 3); 
f(n - 1) = 1/d;
x(n - 1) = f(n - 1) *(lam * y(n - 1) + mu * x(n - 2) - x(n - 3));

mu = 2 - e(n - 2); 
e(n - 1) = mu * f(n - 1);
d = a1 - mu * e(n - 1) - f(n - 2); 
f(n) = 1/d;
x(n) = f(n) *(lam * y(n) + mu * x(n - 1) - x(n - 2));




%Solve the second triangular system and find avg squared error
sq =(y(n) - x(n))^2;
x(n - 1) = x(n - 1) + e(n - 1) * x(n);
sq = sq + (y(n - 1) - x(n - 1))^2;
for j = n - 2 : -1 : 1
  x(j) = x(j) + e(j) * x(j + 1) - f(j) * x(j + 2);
  sq = sq + (y(j) - x(j))^2;
end

sq = sq/n;
%Compute GCV score
g2 = f(n); tr = g2; h = e(n - 1) * g2;
g1 = f(n - 1) + e(n - 1) * h; tr = tr + g1;
for k = n - 2 : -1 : n - nc + 1
  q = e(k) * h - f(k) * g2;
  h = e(k) * g1 - f(k) * h; g2 = g1;
  g1 = f(k) + e(k) * h - f(k) * q;
  tr = tr + g1;
end
tr =(2 * tr - rem(n, 2) * g1) * lam/n;
score = sq/(1 - tr)^2;

end

