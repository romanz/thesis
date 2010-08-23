% [Q] = dirichlet(grid, dir)
function [Q] = dirichlet(grid, dir)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    boundary = find(J); % boundary indices
    Q = sparse(boundary, 1:numel(boundary), ones(size(boundary)), ...
        grid.numel, numel(boundary));
end
