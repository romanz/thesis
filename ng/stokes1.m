% Construct Stokes equation for specified dimension.
function [Gp, Lv, Gv] = stokes1(dim, gridV, gridP)
    Lv = laplacian(gridV.I, gridV.X, gridV.Y);
    % full(Lv)
    dir = (1:numel(gridP.sz)) == dim;
    K = gridV.I | shift(gridV.I, -dir); % add first row/column
    Gv = gradient(K, gridV.X, gridV.Y, dir);
    % full(Gv)
    K = gridP.I & shift(gridP.I, -dir); % remove last row/column
    Gp = gradient(K, gridP.X, gridP.Y, dir);
    % full(Gp)
end
