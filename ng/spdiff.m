% Sparse matrix gradient representation for finite difference operator.
function D = spdiff(grid, dim)
    dir = (1:2 == dim); assert(any(dir));
    J1 = grid.I | shift(grid.I, dir);
    J0 = shift(J1, -dir);
    J = [find(J0) find(J1)]; % column indices
    N = size(J, 1); % # of equations
    I = repmat(1:N, 1, 2); % row indices
    
    % Finite difference operator
    D = sparse(I, J, repmat([-1 1], N, 1), N, grid.numel); 
end
