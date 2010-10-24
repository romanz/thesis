function [W, gridW] = vorticity(Vx, Vy, gridVx, gridVy)

gridW = init_grid(gridVx.x, gridVy.y);
W = (diff(gridVy.X .* Vy, [], 1) ./ diff(gridVy.X, [], 1) - ...
     diff(Vx, [], 2) ./ diff(gridVx.Y, [], 2));
W = W ./ gridW.X;

end