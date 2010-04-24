function [A, f] = subst(A, J, W, v)
% A(k, :) x(:) = 0
% W(k, :) x(J(k, :)) = v(k)
[N, M] = size(A); % M = # of variables, N = # of equations.

S = spdiag(1./ W(:, 1)); % Rescaling according to first coordinate
W = S * W;
v = S * v;

W = [W, v];

P = J(:, 1); J(:, 1) = 0;
I = repmat(P, [1 size(W, 2)]);
J = [J, repmat((M+1), size(J, 1), 1)];
K = logical(J);

Q = (1:M)'; 
Q(P) = [];
I = I(K); 
J = J(K); 
W = W(K);
T = sparse([I(:); Q], [J(:); Q], [W(:); ones(size(Q))], M, M+1);
A = A * T;

f = A(:, end);
A = A(:, 1:end-1);