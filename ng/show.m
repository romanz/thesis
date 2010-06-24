function show(grid, M, msg)
% M(1, 1) = NaN;
% M(1, end) = NaN;
% M(end, 1) = NaN;
% M(end, end) = NaN;
mesh(grid.x, grid.y, M.')
title(msg)
xlabel('R')
ylabel('\theta')
