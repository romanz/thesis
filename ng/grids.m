function [center, interior, xstag, ystag] = grids(x, y)
    x = x(:); y = y(:); 
    % extended grid (for ghost points)
    xg = [2*x(1) - x(2); x; 2*x(end) - x(end-1)];
    yg = [2*y(1) - y(2); y; 2*y(end) - y(end-1)];
    % cell-centered coordinates
    xc = average(xg, [1; 1]/2);
    yc = average(yg, [1; 1]/2);
    
    % NDGRID convention
    center = grid(xc, yc);    
    xstag = grid(xg(2:end-1), yc);
    ystag = grid(xc, yg(2:end-1));    
    interior = grid(xc(2:end-1), yc(2:end-1));
    interior.I = true(size(interior.I)); % no boundary
end
