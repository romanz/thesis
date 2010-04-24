% Substitude W*x(J)=u into Ax=0, giving Bx=f
function [f, A] = subst(A, J, W, u)
% A(k, :) x(:) = 0
% W(k, :) x(J(k, :)) = u(k)
[N, M] = size(A); % M = # of variables, N = # of equations.

S = spdiag(1./ W(:, 1)); % Rescaling according to first coordinate
v = S * v;

P = J(:, 1); J(:, 1) = 0; % P is to-be-eliminated variable vector
f = A * sparse(P, ones(size(P)), -u, M, 1);

if nargout >= 2 % Lazy evaluation for A:
    W = S * W;
    I = repmat(P, [1 size(W, 2)]);
    K = logical(J);

    Q = (1:M)'; Q(P) = []; % Q is P's complement    

    I = I(K); J = J(K); W = W(K);
    T = sparse([I(:); Q], [J(:); Q], [W(:); ones(size(Q))], M, M);
    A = A * T;
end
