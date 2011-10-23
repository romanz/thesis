function asymp

    syms pi b r t a g U1 U2 real
    syms A1 A2 A3 real
    integral = @(f) int(f, 0, pi);
    
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi + log(c); c*Dr(phi) - Dr(c)]);
    slip = @(phi) bnd(dukhin_derjaguin(phi, g));

    phi1 = (1/4 * r^(-2) - r) * cos(t);
    phi2 = (3/32*r^-4  - 3/8*r^-1)*sin(t)^2 - 3/32*r^-4;
    phi3 = (15*U1*a + 6)/64*cos(t) + (-3*U1*a/32)*cos(t)^3 ...
         + (15*U1*a + 6)/64*cos(t)^3/r^2 ...
         + ((12 + 15*U1*a)/64*cos(t) - (27 + 18*U1*a)/96*cos(t)^3)/r^3 ...
         + ((2*U1*a - 1)/64*cos(t) + (3 - 6*U1*a)/64*cos(t)^3)/r^5 ...
         + (27 + 9*U1*a)/576*cos(t)^3/r^6;
    
    c1 = 3/4 * r^(-2) * cos(t);
    c2 = a*U1*3/8*((r^-1 + 1/2*r^-4)*sin(t)^2 - 1/2*r^-4);
    v1 = U1 * [-(1 - r^-3) * cos(t); (1 + (r^-3)/2) * sin(t)];
    v2 = U2 * curl( (cos(t)*sin(t)^2*(r^-2 - 1)) );
    
    phi2 = phi2 + 3*(1/16 - a*U1/32)/r + (a*U1/32 - 1/16)*(3*cos(t)^2 - 1)/r^3;
    c2 = c2 + 3*(1/16 - a*U1/32)/r + ((5*U1*a)/32 + 1/16)*(3*cos(t)^2 - 1)/r^3;
    
    phi = b * phi1 + b^2 * phi2 + b^3 * phi3;
    c = 1 + b * c1 + b^2 * c2;
    v = b * v1 + b^2 * v2;
    p = b^2 * U2 * 2 * r^-3 * (1 - 3*cos(t)^2);
    
    eq1 = divergence(c * gradient(phi));
    assert_zero( series(eq1, b, 0, 3), 'Poisson' )
    
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
    u1 = integral(vs        * sin(t)^2/b) / int(sin(t), 0, pi);
    w1 = integral(bnd(v(2)) * sin(t)^2/b) / integral(sin(t));
    u2 = integral(vs        * b^-2*2*sin(2*t)/pi);
    w2 = integral(bnd(v(2)) * b^-2*2*sin(2*t)/pi);
    save asymp
end

function assert_zero(e, msg)
    e = simple(e);
    z = all(e == 0);
    if ~z
        pretty(e)
        save('assert.mat')
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
