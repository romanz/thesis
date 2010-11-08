function test
    syms a r t U W
    clc
    psi = (U/2) * (r^2 - (3*a*r)/2 + a^3/(2*r))*sin(t)^2;
    V = curl(psi);
    L = simple(vector_laplacian(V));
    P = -3*U*a*cos(t)/(2*r^2);
    gP = gradient(P);
    assert(all(L == gP))
    F = radial_stress(V, P, a);

    psi = (W/2) * (a*r - a^3/r)*sin(t)^2;
    V = curl(psi);
    subs(V, r, a)
    L = simple(vector_laplacian(V));
    P = W*a*cos(t)/r^2;
    gP = gradient(P);
    assert(all(L == gP))
    F = radial_stress(V, P, a);
    
    psi1 = (U/2) * (r^2 - (3*a*r)/2 + a^3/(2*r))*sin(t)^2;
    psi2 = (W/2) * (a*r - a^3/r)*sin(t)^2;
    psi2 = subs(psi2, W, 3*U/2);
    psi = simple(psi1 + psi2);
    pretty(psi)
    V = simple(curl(psi))
    pretty(V)
    L = simple(vector_laplacian(V))    
    P = 0; % Note that the pressure is constant!
    gP = gradient(P);
    assert(all(L == gP))
    F = radial_stress(V, P, a)
    save results
end

function F = gradient(f)
    syms r t
    F = [diff(f, r); diff(f, t) / r];
end

function d = divergence(F)
    syms r t
    d = diff(r^2 * F(1), r) / (r^2) + diff(sin(t) * F(2), t) / (r*sin(t));
end

function L = scalar_laplacian(f)
    L = divergence(gradient(f));
end

function V = curl(psi)
    syms r t
    V = [diff(psi, t) / r; -diff(psi, r)] / (r * sin(t));
end

function [L] = vector_laplacian(F)
    syms r t
    L = [scalar_laplacian(F(1)) ...
         - 2 * F(1) / (r^2) - 2 * diff(sin(t) * F(2), t) / (r^2 * sin(t)); ...
         scalar_laplacian(F(2)) ...
         - F(2) / (r * sin(t))^2 + 2 * diff(F(1), t) / r^2];
end

function F = radial_stress(V, P, a)
    syms r t
    S = [-P + 2*diff(V(1), r);
        diff(V(2), r) + (diff(V(1), t) - V(2)) / r];
    S = subs(S, r, a);
    S = simple(S)
    f = simple([cos(t), -sin(t)] * S)
    F = int(f * 2*pi*a*sin(t), t, 0, pi);
end