% Create sparse matrix for linear intepolation 
% along specified dimension. 
function [P, target] = sparse_interpolate(grid, dim, target)
    dir = (1:2 == dim); assert(any(dir));
    J1 = grid.I | shift(grid.I, dir);
    J0 = shift(J1, -dir);
    J = [find(J0) find(J1)]; % column indices
    N = size(J, 1); % # of equations
    I = repmat(1:N, 1, 2); % row indices
    
    % Interpolation operator
    if nargin < 3
        h = [1 1]/2;
        weights = repmat(h, N, 1);
        x = grid.x(:);
        y = grid.y(:);
        switch dim
            case 1, x = convn(x, h(:), 'valid'); y = y(2:end-1);
            case 2, y = convn(y, h(:), 'valid'); x = x(2:end-1);
        end
        target = init_grid(x, y);
    else
        error('NIY')
    end 
    P = sparse(I, J, weights, N, grid.numel); 
end
