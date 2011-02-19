% G = INIT_GRID(x, y)
%   Create regular grid object for specified (x,y) coordinates.
function G = init_grid(x, y)
    G.x = x(:);
    G.y = y(:);
    [G.X, G.Y] = ndgrid(G.x, G.y);
    G.sz = [numel(x), numel(y)];
    G.numel = prod(G.sz);
    G.I = interior(G.sz);
end
