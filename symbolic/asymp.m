function asymp

    syms pi b r t a g U1 U2 U3 real
    syms U3a U3b real
    
    order = 3;
    a = sym(0);
    
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi + log(c); c*Dr(phi) - Dr(c)]);
    slip = @(phi) bnd(dukhin_derjaguin(phi, g));

    %% O(beta) : Homogeneous equations
    % Boundary conditions (by spherical harmonics)
    phi1 = (1/4 * r^(-2) - r) * cos(t);
    c1 = 3/4 * r^(-2) * cos(t);
    % Streamfunction
    psi1 = (r^-1 - r^2)*sin(t)^2/2;
    assert_zero(E4(psi1), '')
    v1 = U1 * curl( psi1 ); 
    
    %% O(beta^2) : Homogeneous equations
    % Non-homogeneous solution (Phi & C)
    phi2 = (3/32*r^-4  - 3/8*r^-1)*sin(t)^2 - 3/32*r^-4;
    c2 = a*U1*3/8*((r^-1 + 1/2*r^-4)*sin(t)^2 - 1/2*r^-4);
    % Boundary conditions (by spherical harmonics)
    phi2 = phi2 + 3*(1/16 - a*U1/32)/r + (a*U1/32 - 1/16)*(3*cos(t)^2 - 1)/r^3;
    c2 = c2 + 3*(1/16 - a*U1/32)/r + ((5*U1*a)/32 + 1/16)*(3*cos(t)^2 - 1)/r^3;
    % Streamfunction (homogeneous Stokes, for sin(2t) boundary conditions)
    psi2 = (r^-2 - 1)*(cos(t)*sin(t)^2);
    assert_zero(E4(psi2), '')
    v2 = U2 * curl( psi2 );
    p2 = U2 * 2 * r^-3 * (1 - 3*cos(t)^2);

    %% O(beta^3)
    % Non-homogeneous solution (Phi & C)
    phi3 = (15*U1*a + 6)/64*cos(t) + (-3*U1*a/32)*cos(t)^3 ...
         + (15*U1*a + 6)/64*cos(t)^3/r^2 ...
         + ((12 + 15*U1*a)/64*cos(t) - (27 + 18*U1*a)/96*cos(t)^3)/r^3 ...
         + ((2*U1*a - 1)/64*cos(t) + (3 - 6*U1*a)/64*cos(t)^3)/r^5 ...
         + (27 + 9*U1*a)/576*cos(t)^3/r^6;
    c3 = (3*U1^2*a^2*cos(t)^3)/32 - (3*U1*a*cos(t)*(5*U1*a + 2))/64 ...
        - (3*a*cos(t)*(- 8*a*U1^2*cos(t)^2 + 7*a*U1^2 + 2*U1 + 32*U2*cos(t)^2 - 32*U2))/(128*r^3) ...
        + (a*(cos(t) - 3*cos(t)^3)*(5*a*U1^2 + 2*U1 + 16*U2))/(128*r^5) ...
        + (U1^2*a^2*cos(t)^3)/(32*r^6) - (3*U1*a*cos(t)^3*(5*U1*a + 2))/(64*r^2);
    % Boundary conditions (by spherical harmonics)
    phi3 = phi3 + ...
        ((183*U1^2*a^2)/2560 - (839*U1*a)/2560 - 47/256) * cos(t)/r^2 + ...
        (- (19*U1^2*a^2)/5120 + (97*U1*a)/5120 + (3*U2*a)/320 + 21/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    c3 = c3 + ...
        ((797*U1^2*a^2)/2560 + (419*U1*a)/2560 - (U2*a)/5 + 87/1280) * cos(t)/r^2 + ...
        ((59*U1^2*a^2)/5120 + (103*U1*a)/5120 + (69*U2*a)/320 + 3/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    
    % Streamfunction : non-homogenous solution
    psi3 = (sin(t)^2*(-1848*r^6*sin(t)^2 + 2079*r^3*sin(t)^2 - 924*r^3 + 7*sin(t)^2 - 10))/(19712*r^5);
    p3 = (-cos(t)*(- 1584*r^6*cos(t)^2 + 1584*r^6 + 1188*r^3*cos(t)^2 - 396*r^3 + 45*cos(t)^2 + 17))/(1408*r^8);
    % Homogeneous solution for boundary conditions (slip of U3a*sin(t) + U3b*sin(3t))
    A1 = 829/28160 + U3b*2/5;         A2 = -761/5632 - U3b*2;
    A3 = 35/264 - U3a/3 + U3b*5/3;    A4 = -209/3360 + U3a/3 - U3b*1/15;
    psi3 = psi3 - (A1*(1-5*cos(t)^2)*sin(t)^2/r^3 + A2*sin(t)^4/r + ...
                   A3*sin(t)^2/r + A4*r^2*sin(t)^2);
    p3 = p3 - A2*(6*cos(t) - 10*cos(t)^3)/r^4;
    
    % Stokes flow term
    U3 = (2 - U1*a)/16; % zero slip, U3 flow at infinity
    psi3 = psi3 - U3*(r^2 - (3/2)*r + (1/2)/r)*sin(t)^2/2;
    p3 = p3 + U3*3/2*cos(t)/r^2;
    
    v3 = curl( psi3 );
    
    %% Variables
    phi = b * phi1 + b^2 * phi2 + b^3 * phi3;
    c = 1 + b * c1 + b^2 * c2 + b^3 * c3;
    v = b * v1 + b^2 * v2 + b^3 * v3;
    p = b^2 * p2 + b^3 * p3;
    
    %% Equations
    eq1 = divergence(c * gradient(phi));
    assert_zero( series(eq1, b, 0, order), 'ion flux' )
    
    eq2 = scalar_laplacian(c) - a * sum(v .* gradient(c));
    assert_zero( series(eq2, b, 0, order), 'salt flux' )
    
    eq3 = [vector_laplacian(v) - gradient(p) + scalar_laplacian(phi)*gradient(phi); ...
           divergence(v)];

    assert_zero( series(eq3, b, 0, order), 'Stokes' )
    
    eq4 = series( cond(phi, c), b, 0, order );
    assert_zero(eq4, 'coupled boundary condition on R=1 for C and Phi' )
    
    f = force(v, p, phi);
    f = series(f, b, 0, order);
    
    assert_zero( f, 'total force' )
        
    vb = simple(subs(v, r, 1)); % slipping
    
    W1 = 2*log((1+1/sqrt(g))/2);
    W2 = 9/(16*(sqrt(g)+1)) - W1*(W1*a + 1)*3/16;
    
    W3a = 69/(512*(g^(1/2) + 1)) - (123*log((g^(1/2) + 1)/(2*g^(1/2))))/1280 ...
        - 27/(512*(g^(1/2) + 1)^2) - (1407*W1^2*a^2*log((g^(1/2) + 1)/(2*g^(1/2))))/5120 ...
        + (21*W1*a)/(256*(g^(1/2) + 1)) - (1899*W1*a*log((g^(1/2) + 1)/(2*g^(1/2))))/5120 ...
        - (9*W2*a*log((g^(1/2) + 1)/(2*g^(1/2))))/320;
    W3b = (9*log((g^(1/2) + 1)/(2*g^(1/2))))/256 - 27/(512*(g^(1/2) + 1)) ...
        - 27/(512*(g^(1/2) + 1)^2) + (57*W1^2*a^2*log((g^(1/2) + 1)/(2*g^(1/2))))/1024 ...
        - (27*W1*a)/(256*(g^(1/2) + 1)) + (93*W1*a*log((g^(1/2) + 1)/(2*g^(1/2))))/1024 ...
        - (9*W2*a*log((g^(1/2) + 1)/(2*g^(1/2))))/64;
    vs = slip(phi);
    vs = series(vs, b, 0, order);
    slip_error = subs(vs - vb(2), [U1 U2 U3a U3b], [W1 W2 W3a W3b]);
    
    assert_zero( vb(1), 'no penetration' );
    assert_zero( slip_error, 'slip condition' );

    assert_zero(limit(gradient(phi), r, inf) - b*[-cos(t);sin(t)], 'boundary R=Inf for Phi' )
    assert_zero(limit(c, r, inf) - 1, 'boundary R=Inf for C' ) % ??? %
    
    Uinf = simple([-cos(t), sin(t)]*simple(limit(v, r, inf)));
    Uinf = subs(Uinf, [U1 U2 U3a U3b], [W1 W2 W3a W3b]);
    
    if a == 0
        Winf = b*(1 - 11/320*b^2)*W1 ...
             - (9/320*(sqrt(g)+1)^-2 - 31/320*(sqrt(g)+1)^-1 - 1/1680)*b^3;
        assert_zero(Uinf - Winf, 'Uinf (for a=0)');
    end
    
    save asymp
    
    fprintf('Done.\n')
end

function assert_zero(e, msg)
    e = simple(e);
    z = all(e == 0);
    if msg
        fprintf('Checking %s... ', msg);
    end
    if ~z
        pretty(e)
        save('assert.mat')
        error('assert:zero', msg)        
    end
    if msg
        fprintf('OK.\n');
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
    syms r t b pi
    Vr = v(1);
    Vt = v(2);
    Fr = simple( -p + 2*Dr(Vr) + ((Dr(phi)^2 - (Dt(phi)/r)^2))/2 );
    Ft = simple( Dr(Vt) + (Dt(Vr) - Vt)/r + Dr(phi)*Dt(phi)/r    );
    df = Fr*cos(t) - Ft*sin(t);
    df = subs(df, r, 1);
    df = simple(df);
    f = int(df * 2 * pi * sin(t), t, 0, pi);
end
