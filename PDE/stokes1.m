% Construct Stokes equation for specified dimension.
function [Gp, Lv, Gv] = stokes1(sz, X, Y, dim)
    [Xs, Ys, dir] = staggered_grid(sz, X, Y, dim); 
    I = interior(sz - dir); 
    Lv = laplacian(I, Xs, Ys);
    Lv = Lv(I, :);
    % full(Lv)

    K = I | shift(I, -dir); % add first row/column
    Gv = gradient(K, Xs, Ys, dir);
    % full(Gv)

    I = interior(sz);
    K = I & shift(I, -dir); % remove last row/column
    Gp = gradient(shave(K, 1, 1), X(I), Y(I), dir);
    % full(Gp)
end
