function V = curl(psi)
    syms r t
    V = [diff(psi, t) / r; -diff(psi, r)] / (r * sin(t));
end
