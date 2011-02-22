function S = stokes(grid)

    [Gp_x, Lv_x0, Dv_x] = stokes1d(1, grid.Vx, grid.P);
    [Gp_y, Lv_y0, Dv_y] = stokes1d(2, grid.Vy, grid.P);

    Ix = grid.Vx.I;
    Iy = grid.Vy.I;

    % -2Vr / r^2
    Lv_x1 = expand(grid.Vx.I)' * spdiag(-2./grid.Vx.X.^2);

    % -2/(r^2 sin(theta)) d(Vtheta sin(theta))/dtheta
    X = average(grid.Vy.X, [1 1]/2);
    Xi = grid.Vx.X(:, 2:end-1);
    J = grid.Vx.I(:, 2:end-1);
    Lv_x2 = spdiag(-2./( sin(grid.Vx.Y(Ix)) .* grid.Vx.X(Ix).^2 )) * ...
            expand( J )' * interpolator(X, Xi) * ...
            spderiv(grid.Vy.Y, 2) * spdiag( sin(grid.Vy.Y) );

    % -Vtheta/(r * sin theta)^2
    Lv_y1 = expand(grid.Vy.I)' * spdiag(-(grid.Vy.X .* sin(grid.Vy.Y) ).^(-2));

    % 2/r^2 dVr/dtheta
    Y = average(grid.Vx.X, [1 1]/2);
    Yi = grid.Vy.X(2:end-1, :);
    J = grid.Vy.I(2:end-1, :);
    Lv_y2 = spdiag( 2./grid.Vy.X(Iy).^2 ) * expand( J )' * ...
            interpolator(Y, Yi) * spderiv(grid.Vx.Y, 2);

    Z = sparse(grid.P.numel, grid.P.numel);

    S = [-Lv_x0 - Lv_x1, -Lv_x2,         Gp_x; ...
         -Lv_y2,         -Lv_y0 - Lv_y1, Gp_y; ...
         Dv_x,            Dv_y,          Z];

end

% Construct Stokes equation for specified dimension.
function [Gp, Lv, Dv] = stokes1d(dim, gridV, gridP)
    [D1, G1] = operators(gridV, 1);
    [D2, G2] = operators(gridV, 2);    
    Lv = D1 * G1 + D2 * G2;
    
    Dv = sparse_divergence(gridV, dim, gridP);    
    Gp = sparse_gradient(gridP, dim);
end

function D = spderiv(G, dim)
    D = spdiff(true(size(G)), dim);
    D = spdiag(1 ./ (D * G(:))) * D;
end
