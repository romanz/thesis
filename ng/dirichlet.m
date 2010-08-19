% [P, Q] = dirichlet(grid, dir)
function [P, Q] = dirichlet(grid, dir)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    boundary = find(J); % boundary indices
    P = sparse(grid.numel, grid.numel); 
    Q = sparse(boundary, 1:numel(boundary), ones(size(boundary)), ...
        grid.numel, numel(boundary));
end
