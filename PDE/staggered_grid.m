function [Xs, Ys, dir] = staggered_grid(sz, X, Y, dim)
dir = double((1:numel(sz)) == dim);
h = ones(dir + 1) / 2;
Xs = average(X, h);
Ys = average(Y, h);
