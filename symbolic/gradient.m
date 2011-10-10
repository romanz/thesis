function F = gradient(f)
    syms r t
    F = [diff(f, r); diff(f, t) / r];
end
