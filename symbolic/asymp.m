function asymp

    syms pi b r t A g
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi + log(c); c*Dr(phi) - Dr(c)]);
    slip = @(phi) bnd(dukhin_derjaguin(phi, g));

    phi1 = b * (1/4 * r^(-2) - r) * cos(t);
    c1 = 1 + b * 3/4 * r^(-2) * cos(t);
    taylor( cond(phi1, c1), 2, b )

    v1 = slip(phi1);
    v1 = simplify(taylor(v1, 2, b));
    u1 = int(v1 * sin(t)^2, 0, pi) / int(sin(t), 0, pi);
end

function DfDr = Dr(f)
    DfDr = diff(f, 'r');
end

function DfDt = Dt(f)    
    DfDt = diff(f, 't');
end

function v = dukhin_derjaguin(phi, g)
    v = 4 * log((exp((-phi-log(g))/2) + 1)/2) * Dt(phi);
end
