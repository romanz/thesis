function g = E2(f)
    syms r t
    g = diff(f, r, 2) + (sin(t)/r^2)*diff(diff(f, t)/sin(t), t);
end
