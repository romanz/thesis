function [u, residuals] = iterate(u, A, f, R, iters, iter_type)
% Apply an iterative scheme to solve u' = Cu + d.

if nargin < 4
    iter_type = ''; % Use plain Jacobi iteration.
end

% Keep the original size and convert to column vector:
sz = size(u); 
u = u(:);
N = numel(u);
T = speye(N) - R * A;
d = R * f;
residuals = zeros(iters, 1);
h = progress([], 0, sprintf('%s Iteration', iter_type));

if strcmpi(iter_type, 'redblack')
    % Prepare Checkerboard pattern
    P = cumsum(ones(sz), 1) + cumsum(ones(sz), 2);
    P = logical(mod(P(:), 2)); 
    % Create logical index Red and Black matrices
    red = P;
    black = ~P;
    % Split {T,d} into their red and black version:
    T_red = T(red, :);
    d_red = d(red, :);    
    T_black = T(black, :);
    d_black = d(black, :);
    % Iterate using Red/Black Gauss-Seidel method:
    for iter = 1:iters
        u_old = u;
        u(red) = T_red * u + d_red; % Update u's red interior.
        u(black) = T_black * u + d_black; % Update u's black interior.
        residuals(iter) = norm(u - u_old);
        h = progress(h, iter/iters);
    end
else
    % Plain Jacobi iteration:
    for iter = 1:iters
        u_old = u;
        u = T * u + d;
        residuals(iter) = norm(u - u_old);
        h = progress(h, iter/iters);
    end
end
progress(h, []);
u = reshape(u, sz);
