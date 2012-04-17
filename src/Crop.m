function op = Crop(op, d)
    g = op.grid;
    r = g.r(1+d(1):end-d(1));
    t = g.t(1+d(2):end-d(2));
    g = Grid(r, t);
    op = Selector(g, op);
end
