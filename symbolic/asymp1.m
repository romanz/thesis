function asymp

    syms pi b r t a g U1 U2 U3 real
    syms U3a U3b real

    a = 0;
    order = 3;
        
    bnd = @(f) subs(f, r, 1);
    cond = @(phi, c) bnd([phi + log(c); c*diff(phi, r) - diff(c, r)]);

    %% O(beta) : Homogeneous equations
    % Boundary conditions (by spherical harmonics)
    phi1 = (1/4 * r^(-2) - r) * cos(t);
    c1 = 3/4 * r^(-2) * cos(t);
    
    %% O(beta^2) : Homogeneous equations
    % Non-homogeneous solution (Phi & C)
    phi2 = (3/32*r^-4  - 3/8*r^-1)*sin(t)^2 - 3/32*r^-4;
    % Boundary conditions (by spherical harmonics)
    phi2 = phi2 + 3*(1/16)/r + (-1/16)*(3*cos(t)^2 - 1)/r^3;
    c2 = 3*(1/16)/r + (1/16)*(3*cos(t)^2 - 1)/r^3;

    %% O(beta^3)
    % Non-homogeneous solution (Phi & C)
    phi3 = (6)/64*cos(t) ...
         + (6)/64*cos(t)^3/r^2 ...
         + ((12)/64*cos(t) - (27)/96*cos(t)^3)/r^3 ...
         + ((-1)/64*cos(t) + (3)/64*cos(t)^3)/r^5 ...
         + (27)/576*cos(t)^3/r^6;
    % Boundary conditions (by spherical harmonics)
    phi3 = phi3 + ...
        (- 47/256) * cos(t)/r^2 + ...
        (21/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    c3 = (87/1280) * cos(t)/r^2 + ...
         (3/1280) * (5*cos(t)^3 - 3*cos(t))/r^4;
    
    %% Variables
    phi = b * phi1 + b^2 * phi2 + b^3 * phi3;
    c = 1 + b * c1 + b^2 * c2 + b^3 * c3;
    
    %% Equations
    eq1 = divergence(c * grad(phi));
    assert_zero( series(eq1, b, 0, order), 'ion flux' )
    
    eq2 = scalar_laplacian(c);
    assert_zero( series(eq2, b, 0, order), 'salt flux' )
    
    eq4 = series( cond(phi, c), b, 0, order );
    assert_zero(eq4, 'coupled boundary condition on R=1 for C and Phi' )
    
    assert_zero(limit(grad(phi), r, inf) - b*[-cos(t);sin(t)], 'boundary R=Inf for Phi' )
    
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
