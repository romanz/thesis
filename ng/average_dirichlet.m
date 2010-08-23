function [P, Q] = average_dirichlet(grid, dir)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    I = shift(J, -dir); % interior neighbourhood logical mask
    boundary = find(J); % boundary indices
    interior = find(I); % interior indices
    P = sparse(boundary, interior, repmat(-1, size(boundary)), ...
        grid.numel, grid.numel) + speye(grid.numel);
    Q = sparse(boundary, 1:numel(boundary), repmat(2, size(boundary)), ...
        grid.numel, numel(boundary));
end
