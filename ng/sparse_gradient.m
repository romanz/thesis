% Create sparse matrix gradient representation G, 
% computed on interior cells' edges.
function D = sparse_gradient(grid, dim)
    dir = (1:2 == dim);
    J1 = grid.I | shift(grid.I, dir);
    J0 = shift(J1, -dir);
    J = [find(J0) find(J1)]; % column indices
    N = size(J, 1); % # of equations
    I = repmat(1:N, 1, 2); % row indices
    
    % Finite difference operator
    D = sparse(I, J, repmat([-1 1], N, 1), N, grid.numel); 
    switch dim
        case 1, g = grid.X; % ds = d{r}
        case 2, g = grid.Y .* grid.X; % ds = r d{theta}
        otherwise, error('invalid dimension');
    end
    dg = D * g(:); % Scale by grid differences reciprocal.
    D = spdiag(1./dg) * D;
end

