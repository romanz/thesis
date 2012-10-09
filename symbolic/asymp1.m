function asymp1

    syms pi b r t a g U1 U2 U3 real
    a = 0;
    order = 3;
        
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi - log(c); diff(phi + log(c), r)]);

    %% O(beta) : Homogeneous equations
    % Boundary conditions (by spherical harmonics)
    phi1 = 1 - (1/2) * r^(-1);
    c1 = (1/2) * r^(-1);
    
    %% O(beta^2) : Homogeneous equations
    % Boundary conditions (by spherical harmonics)
    a1 = 1/16; a3 = 3/16;
%     syms a1 a2 a3 a4 real
    phi2 = a1/r;% + a2*(3*cos(t)^2 - 1)/r^3;
    c2 = a3/r;% + a4*(3*cos(t)^2 - 1)/r^3;

    %% O(beta^3)
    % Boundary conditions (by spherical harmonics)
    a1 = 1/192; a3 = 11/192; 
%     syms a1 a2 a3 a4 real
    phi3 = a1/r; % * cos(t)/r^2 + a2 * (5*cos(t)^3 - 3*cos(t))/r^4;
    c3 = a3/r; % * cos(t)/r^2 + a4 * (5*cos(t)^3 - 3*cos(t))/r^4;
    
    %% Variables
    phi = b * phi1 + b^2 * phi2 + b^3 * phi3;
    c = 1 + b * c1 + b^2 * c2 + b^3 * c3;
    
    %% Equations
    eq1 = divergence(grad(phi));
    assert_zero( series(eq1, b, 0, order), 'ion flux' )
    
    eq2 = scalar_laplacian(c);
    assert_zero( series(eq2, b, 0, order), 'salt flux' )
    
    eq4 = series( cond(phi, c), b, 0, order );
    assert_zero(eq4, 'coupled boundary condition on R=1 for C and Phi' )
    
    assert_zero(limit(phi, r, inf) - b, 'boundary R=Inf for Phi' )
    
    if a == 0
        assert_zero(limit(c, r, inf) - 1, 'boundary R=Inf for C' ) 
    end
    
    save asymp1
    
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

% inline(vectorize(char
function V = curl(psi)
    syms r t
    V = [diff(psi, t) / (r^2 * sin(t)); ...
        -diff(psi, r) / (r * sin(t))];
    V = simple(V);
end

function [L] = vector_laplacian(F)
    syms r t
    L = [scalar_laplacian(F(1)) ...
         - 2 * F(1) / (r^2) - 2 * diff(sin(t) * F(2), t) / (r^2 * sin(t)); ...
         scalar_laplacian(F(2)) ...
         - F(2) / (r * sin(t))^2 + 2 * diff(F(1), t) / r^2];
end

function L = scalar_laplacian(f)
    L = simplify(divergence(grad(f)));
end

function F = grad(f)
    syms r t
    F = [diff(f, r); diff(f, t) / r];
end

function d = divergence(F)
    syms r t
    d = diff(r^2 * F(1), r) / (r^2) + diff(sin(t) * F(2), t) / (r*sin(t));
end
