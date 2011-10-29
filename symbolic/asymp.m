function asymp

    syms pi b r t a g U1 U2 U3 real
    syms A1 A2 A3 A4 real
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
    c3 = (3*U1^2*a^2*cos(t)^3)/32 ...
        - (3*a*cos(t)*(- 8*a*U1^2*cos(t)^2 + 7*a*U1^2 + 2*U1 + 32*U2*cos(t)^2 - 32*U2))/(128*r^3) ...
        - (3*U1*a*cos(t)*(5*U1*a + 2))/64 ...
        + (a*(cos(t) - 3*cos(t)^3)*(5*a*U1^2 + 2*U1 + 16*U2))/(128*r^5) ...
        + (U1^2*a^2*cos(t)^3)/(32*r^6) - (3*U1*a*cos(t)^3*(5*U1*a + 2))/(64*r^2);
    
    v1 = U1 * curl( (r^-1 - r^2)*sin(t)^2/2 ); 
    v2 = U2 * curl( (r^-2 - 1)*(cos(t)*sin(t)^2) );
    v3 = curl( - (sin(t)^2*(1848*r^6*sin(t)^2 - 2079*r^3*sin(t)^2 + 924*r^3 - 7*sin(t)^2 + 10))/(19712*r^5));
%     v3 = v3 - curl( (sin(t)^2*(-934 + 238*sin(t)^2))/(19712*r) );
    v3 = v3 - curl(U3*(r^2 - (3/2)*r + (1/2)/r)*sin(t)^2/2);
    
    p2 = U2 * 2 * r^-3 * (1 - 3*cos(t)^2);
    p3 = (-cos(t)*(- 1584*r^6*cos(t)^2 + 1584*r^6 + 1188*r^3*cos(t)^2 - 396*r^3 + 45*cos(t)^2 + 17))/(1408*r^8);
%     p3 = p3 - 17/704*(3*cos(t) - 5*cos(t)^3)/r^4;
    p3 = p3 + U3*3/2*cos(t)/r^2;
    
    phi2 = phi2 + 3*(1/16 - a*U1/32)/r + (a*U1/32 - 1/16)*(3*cos(t)^2 - 1)/r^3;
    c2 = c2 + 3*(1/16 - a*U1/32)/r + ((5*U1*a)/32 + 1/16)*(3*cos(t)^2 - 1)/r^3;
    
    phi3 = phi3 + ...
        ((183*U1^2*a^2)/2560 - (839*U1*a)/2560 - 47/256) * cos(t)/r^2 + ...
        (- (19*U1^2*a^2)/5120 + (97*U1*a)/5120 + (3*U2*a)/320 + 21/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    c3 = c3 + ...
        ((797*U1^2*a^2)/2560 + (419*U1*a)/2560 - (U2*a)/5 + 87/1280) * cos(t)/r^2 + ...
        ((59*U1^2*a^2)/5120 + (103*U1*a)/5120 + (69*U2*a)/320 + 3/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    
    phi = b * phi1 + b^2 * phi2 + b^3 * phi3;
    c = 1 + b * c1 + b^2 * c2 + b^3 * c3;
    v = b * v1 + b^2 * v2 + b^3 * v3;
    p = b^2 * p2 + b^3 * p3;
    
    eq1 = divergence(c * gradient(phi));
    assert_zero( series(eq1, b, 0, 3), 'Poisson' )
    
    eq2 = scalar_laplacian(c) - a * sum(v .* gradient(c));
    assert_zero( series(eq2, b, 0, 3), 'Advection' )
    
    eq3 = [vector_laplacian(v) - gradient(p) + ...
            scalar_laplacian(phi)*gradient(phi); ...
           divergence(v)];

    assert_zero( series(eq3, b, 0, 3), 'Stokes' )
    
    eq4 = series( cond(phi, c), b, 0, 3 );
    assert_zero(eq4, 'Boundary Phi/C' )
    
    f = force(v, p, phi);
    f = series(f, b, 0, 3);
    
    vb = simple(subs(v, r, 1));
    assert_zero( vb(1), 'No penetration' );
    
    W1 = 2*log((1+1/sqrt(g))/2);
    W2 = 9/(16*(sqrt(g)+1)) - W1*(W1*a + 1)*3/16;
    vs = slip(phi);
    vs = series(subs(vs, [U1 U2], [W1 W2]), b, 0, 3);
    save asymp
    fprintf('All OK!\n')
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

function DfDr = Dr(f)
    DfDr = diff(f, 'r');
end

function DfDt = Dt(f)    
    DfDt = diff(f, 't');
end

function v = dukhin_derjaguin(phi, g)
    v = 4 * log((exp((-phi-log(g))/2) + 1)/2) * Dt(phi);
end

function f = force(v, p, phi)
    Vr = v(1);
    Vt = v(2);
    syms r t pi
    Fr = -p + 2*Dr(Vr) + ((Dr(phi)^2 - (Dt(phi)/r)^2))/2;
    Ft = Dr(Vt) + (Dt(Vr) - Vt)/r + Dr(phi)*Dt(phi)/r;
    df = Fr*cos(t) - Ft*sin(t);
    df = subs(df, r, 1);
    df = simple(df);
    f = int(df * 2 * pi * sin(t), t, 0, pi);
end