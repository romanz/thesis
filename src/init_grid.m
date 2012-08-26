% G = INIT_GRID(x, y)
%   Create regular grid object for specified (x,y) coordinates.
function G = init_grid(r, t, boundary)
    if nargin < 3
        boundary = 1;
    end
    G.r = r(:);
    G.t = t(:);
    [G.R, G.T] = ndgrid(G.r, G.t);
    G.sz = [numel(r), numel(t)];
    G.numel = prod(G.sz);
    G.I = true(G.sz);
    if boundary
        G.I([1 end], :) = false;
        G.I(:, [1 end]) = false;
    end
    G.boundary = boundary;
    G.mesh = @(v) mesh(G.R, G.T, v);
end
