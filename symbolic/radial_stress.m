function F = radial_stress(V, P, a)
    syms r t
    S = [-P + 2*diff(V(1), r);
        diff(V(2), r) + (diff(V(1), t) - V(2)) / r];
    S = subs(S, r, a);
    S = simple(S)
    f = simple([cos(t), -sin(t)] * S)
    F = int(f * 2*pi*a*sin(t), t, 0, pi);
end
