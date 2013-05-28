function peclet_test

%% Analytical
% (-y' + u y)' = 0
% y(0) = 0
% y(1) = 1

syms x u
y = (exp(u*x) - 1) ./ (exp(u) - 1);

flux = -diff(y, x) + u*y;
assert( diff(flux, x) == 0 )
assert( subs(y, x, 0) == 0 )
assert( subs(y, x, 1) == 1 )

N = 1001;
X = linspace(0, 1, N);
U = 20;
Y = subs(y, {x, u}, {X, U});

%% Numerical
h = 1/(N-1);
L = (sparse(1:N-2, 1:N-2, +1, N-2, N) + ...
     sparse(1:N-2, 3:N  , +1, N-2, N) + ...
     sparse(1:N-2, 2:N-1, -2, N-2, N));
 
D = (sparse(1:N-2, 1:N-2, -1, N-2, N) + ...
     sparse(1:N-2, 2:N-1  , +1, N-2, N));

eqn = -L/h^2 + U*D/(h);
bnd = sparse([1 2], [1 N], [1 1], 2, N);

A = [eqn; bnd];
b = sparse(N, 1, +1, N, 1);

sol = A \ b;

%% Visual
subplot 211
plot(X, Y, '-', X, sol, '.')

subplot 212
plot(X, Y(:) - sol(:), 'rs')

end
