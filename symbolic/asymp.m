function asymp

    syms pi b r t a g U
    integral = @(f) int(f, 0, pi);
    
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi + log(c); c*Dr(phi) - Dr(c)]);
    slip = @(phi) bnd(dukhin_derjaguin(phi, g));

    phi1 = (1/4 * r^(-2) - r) * cos(t);
    phi2 = (3/32*r^-4  - 3/8*r^-1)*sin(t)^2 - 3/32*r^-4;
    c1 = 3/4 * r^(-2) * cos(t);
    c2 = a*U*3/8*((r^-1 + 1/2*r^-4)*sin(t)^2 - 1/2*r^-4);
    v1 = U * [-(1 - r^-3) * cos(t); (1 + (r^-3)/2) * sin(t)];
    v2 = 0;
    
    phi2 = phi2 + 3*(1/16 - a*U/32)/r + (a*U/32 - 1/16)*(3*cos(t)^2 - 1)/r^3;
    c2 = c2 + 3*(1/16 - a*U/32)/r + ((5*U*a)/32 + 1/16)*(3*cos(t)^2 - 1)/r^3;
    
    phi = b * phi1 + b^2 * phi2;
    c = 1 + b * c1 + b^2 * c2;
    v = b * v1;
    p = 0;
    
    eq1 = divergence(c * gradient(phi));
    assert_zero( series(eq1, b, 0, 2), 'Poisson' )
    
    eq2 = scalar_laplacian(c) - a * sum(v .* gradient(c));
    assert_zero( series(eq2, b, 0, 2), 'Advection' )
    
    eq3 = [vector_laplacian(v) - gradient(p) + ...
            scalar_laplacian(phi)*gradient(phi); ...
           divergence(v)];
            
    assert_zero( series(eq3, b, 0, 2), 'Stokes' )
    
    eq4 = series( cond(phi, c), b, 0, 2 );
    assert_zero(eq4, 'Boundary Phi/C' )

    vs = slip(phi);
    vs = series(vs, b, 0, 2);
    us = int(vs * sin(t)^2, 0, pi) / int(sin(t), 0, pi);
    Us = integral(bnd(v(2)) * sin(t)^2) / integral(sin(t));
    save asymp
end

function assert_zero(v, msg)
    v = simplify(v);
    z = all(v == 0);
    if ~z
        pretty(v)
        error('assert:zero', msg)
    end
end

function f = series(f, x, x0, n)
    f = simplify(taylor(f, n+1, x, x0));
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
