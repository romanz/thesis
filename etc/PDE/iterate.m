function [u, residuals] = iterate(u, T, d, iters, iter_type, msg)
% Apply an iterative scheme to solve u' = Cu + d.

if nargin < 5
    iter_type = ''; % Use plain Jacobi iteration.
end

if nargin < 6
    msg = 'progress';
end

% Keep the original size and convert to column vector:
sz = size(u); 
u = u(:);
N = numel(u);
residuals = zeros(iters, 1);

if strcmpi(msg, 'progress')
    h = progress([], 0, sprintf('%s Iteration', iter_type));
end

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
        if strcmpi(msg, 'progress')
            h = progress(h, iter/iters);
        end
    end
else
    % Plain Jacobi iteration:
    for iter = 1:iters
        u_old = u;
        u = T * u + d;
        residuals(iter) = norm(u - u_old);
        if strcmpi(msg, 'progress')
            h = progress(h, iter/iters);
        end
    end
end
if strcmpi(msg, 'progress')
    progress(h, []);
end
u = reshape(u, sz);
