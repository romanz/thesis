% G = INIT_GRID(x, y)
%   Create regular grid object for specified (x,y) coordinates.
function G = init_grid(x, y, boundary)
    if nargin < 3
        boundary = 1;
    end
    G.x = x(:);
    G.y = y(:);
    [G.X, G.Y] = ndgrid(G.x, G.y);
    G.sz = [numel(x), numel(y)];
    G.numel = prod(G.sz);
    G.I = true(G.sz);
    if boundary
        G.I([1 end], :) = false;
        G.I(:, [1 end]) = false;
    end
    G.boundary = boundary;
end
