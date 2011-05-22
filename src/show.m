function show(func, sol, name, Rmax, varargin)
    g = sol.grid.(name);
    Z = (sol.(name));
    I = g.x < Rmax;
    X = g.X .* cos(g.Y);
    Y = g.X .* sin(g.Y);
    func(X(I, :), Y(I, :), Z(I, :), varargin{:});
end
