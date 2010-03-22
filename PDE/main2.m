function main2
% randn('state', 0);
% fprintf('\n'); 
sz = 1+2.^([1 -inf]*5);

%% Create the grid
x = linspace(0, 2, sz(1));
y = linspace(-1, 1, sz(2));
[X, Y] = ndgrid(x, y);

%% Solve the coupled problem
U = 0 + 0*(rand(sz));
for t = 1:100
[Ac, fc] = C_prob(sz, X, Y, U(2), 1);
C = [U(2); Ac \ fc; 1];
[Au, fu] = U_prob(sz, X, Y, C, C(2), 1);
U = [C(2); Au \ fu; 1];
end
plot(x, C, 'b.-', x, U, 'r-.'); legend('C', 'U'); axis([0 1 0 1])

function [A, f] = C_prob(sz, X, Y, Cl, Cr)
[L, M, I] = laplacian(sz, X, Y);
A = dinv(M) * L;
f = zeros(sz); % Right-Hand Side
f(~I) = NaN; % Not defined on the boundary

Bl = boundary(I, [-1 0]);
Br = boundary(I, [+1 0]);
Bd = boundary(I, [0 -1]);
Bu = boundary(I, [0 +1]);
[A, f] = dirichlet(A, f, find(Bl), Cl);
[A, f] = dirichlet(A, f, find(Br), Cr);
[A, f] = dirichlet(A, f, find(Bd), zeros(nnz(Bd), 1));
[A, f] = dirichlet(A, f, find(Bu), zeros(nnz(Bd), 1));

% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(all(isnan(f(~I))) == 1);
assert(nnz(isnan(A(I, I))) == 0);
assert(nnz(isnan(f(I))) == 0);
% Restrict the operator and RHS on the interior variables
A = A(I, I); 
f = f(I);

function [A, f] = U_prob(sz, X, Y, C, Ul, Ur)
[L, M, I] = laplacian(sz, X, Y, C);
A = dinv(M) * L;
f = zeros(sz); % Right-Hand Side
f(~I) = NaN; % Not defined on the boundary

Bl = boundary(I, [-1 0]);
Br = boundary(I, [+1 0]);
Bd = boundary(I, [0 -1]);
Bu = boundary(I, [0 +1]);
[A, f] = dirichlet(A, f, find(Bl), Ul);
[A, f] = dirichlet(A, f, find(Br), Ur);
% [A, f] = dirichlet(A, f, find(Bd), U0(Bd));
% [A, f] = dirichlet(A, f, find(Bu), U0(Bu));

% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(all(isnan(f(~I))) == 1);
assert(nnz(isnan(A(I, I))) == 0);
assert(nnz(isnan(f(I))) == 0);
% Restrict the operator and RHS on the interior variables
A = A(I, I); 
f = f(I);
