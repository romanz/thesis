% Neumann boundary conditions
% Usage:
% [P, Q] = neumann(grid, dir)
%   A = L * P * expand(gridPhi.I)
%   rhs = -L * Q * vals
function [P, Q] = neumann(grid, dir, extend)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    if nargin > 2 && extend
        d = ~dir;
        J = J | shift(J, +d) | shift(J, -d);
    end
    I = shift(J, -dir); % interior neighbourhood logical mask
    boundary = find(J); % boundary indices
    interior = find(I); % interior indices
    h = [grid.X(J) - grid.X(I), grid.Y(J) - grid.Y(I)];
    h = h * dir(:);
    P = sparse(boundary, interior, ones(size(boundary)), ...
        grid.numel, grid.numel) + speye(grid.numel);
    Q = sparse(boundary, 1:numel(boundary), h, ...
        grid.numel, numel(boundary));
end
