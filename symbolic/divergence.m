function d = divergence(F)
    syms r t
    d = diff(r^2 * F(1), r) / (r^2) + diff(sin(t) * F(2), t) / (r*sin(t));
end
