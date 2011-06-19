function show(func, sol, name, Rmax, dir, varargin)
    g = sol.grid.(name);
    Z = (sol.(name));
    I = g.x < Rmax;
    X = g.X .* cos(g.Y);
    Y = g.X .* sin(g.Y);
    X = X(I, :);
    Y = Y(I, :) * dir;
    Z = Z(I, :);
    func(X, Y, Z, varargin{:});
end
