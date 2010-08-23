function G = grid(x, y)
    G.x = x(:);
    G.y = y(:);
    [G.X, G.Y] = ndgrid(G.x, G.y);
    G.sz = [numel(x), numel(y)];
    G.numel = prod(G.sz);
    G.I = interior(G.sz);
end
