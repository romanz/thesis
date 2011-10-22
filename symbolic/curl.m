function V = curl(psi)
    syms r t
    V = [diff(psi, t) / (r^2 * sin(t)); ...
        -diff(psi, r) / (r * sin(t))];
    V = simple(V);
end
