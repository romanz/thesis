% Construct Stokes equation for specified dimension.
function [Gp, Lv, Dv] = stokes1(dim, gridV, gridP)
    Lv = laplacian(gridV.I, gridV.X, gridV.Y);
    % full(Lv)
    dir = (1:numel(gridP.sz)) == dim;
    K = gridV.I | shift(gridV.I, -dir); % add first row/column
    Dv = polar_divergence(K, gridV.X, gridV.Y, dir);
    % full(Gv)
    K = gridP.I & shift(gridP.I, -dir); % remove last row/column
    Gp = polar_gradient(K, gridP.X, gridP.Y, dir);
    % full(Gp)
end
