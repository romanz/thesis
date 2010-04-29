% Construct Stokes equation for specified dimension.
function [Gp, Lv, Gv] = stokes1(dim, Xs, Ys, X, Y)
    sz = size(0*Xs + 0*Ys);
    I = interior(sz);
    Lv = laplacian(I, Xs, Ys);
    Lv = Lv(I, :);
    % full(Lv)
    dir = (1:numel(sz)) == dim;
    K = I | shift(I, -dir); % add first row/column
    Gv = gradient(K, Xs, Ys, dir);
    % full(Gv)
    I = true(size(0*X+0*Y));
    K = I & shift(I, -dir); % remove last row/column
    Gp = gradient(K, X, Y, dir);
    % full(Gp)
end
