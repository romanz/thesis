function F = grad(f)
    syms r t
    F = [diff(f, r); diff(f, t) / r];
end
