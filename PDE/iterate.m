function [v, residuals] = iterate(v, A, f, R, iters, iter_type)
% Apply an iterative scheme to solve v' = Cv + d.

if nargin < 4
    iter_type = ''; % Use plain Jacobi iteration.
end

% Keep the original size and convert to column vector:
sz = size(v); 
v = v(:);
N = numel(v);
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
        v_old = v;
        v(red) = T_red * v + d_red; % Update v's red interior.
        v(black) = T_black * v + d_black; % Update v's black interior.
        residuals(iter) = norm(v - v_old);
        h = progress(h, iter/iters);
    end
else
    % Plain Jacobi iteration:
    for iter = 1:iters
        v_old = v;
        v = T * v + d;
        residuals(iter) = norm(v - v_old);
        h = progress(h, iter/iters);
    end
end
progress(h, []);
v = reshape(v, sz);
