function [v] = jacobi(A, v, f, interior, T, type)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% J := D^{-1}
% v' = v + J (f - Av) = (I - JA)v + Jf = Cv + d

N = numel(v);
D = diag(A);
f = f(interior); % Restrict RHS.
A = A(interior, :); % Restric operator.

k = nnz(interior);
J = sparse(1:k, 1:k, 1 ./ D(interior)); % Restrict the inverse.
I = speye(N);
% v' = Cv + d
% C := I - JA
% d := Jf
C = sparse(N, N);
d = zeros(N, 1);
C(interior, :) = I(interior, :) - J * A;
d(interior) = J * f(:);
d(~interior) = v(~interior); % Dirichlet boundary conditions

% Keep the original size and convert to column vectors:
sz = size(v); 
v = v(:);

if nargin < 6
    type = ''; % Use plain Jacobi iteration.
end
if strcmpi(type, 'redblack')
    % Prepare Checkerboard pattern
    P = cumsum(ones(sz), 1) + cumsum(ones(sz), 2);
    P = logical(mod(P(:), 2)); 
    % Create logical index Red and Black matrices
    red = P;
    black = ~P;
    % Split {C,d} into their red and black version:
    C_red = C(red, :);
    d_red = d(red);    
    C_black = C(black, :);
    d_black = d(black);
    % Iterate:
    for t = 1:T/2
        v(red) = C_red * v + d_red; % Update v's red interior.
        v(black) = C_black * v + d_black; % Update v's black interior.
    end
else
    for t = 1:T
        v = C * v + d; % Update v's interior.
    end
end
v = reshape(v, sz);
