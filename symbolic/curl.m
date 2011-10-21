function V = curl(psi)
    syms r t
    V = [diff(psi * sin(t), t) / (r * sin(t)); ...
        -diff(psi * r, r)] / r;
end
