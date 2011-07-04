function S = interpolate_matrix(x, xi)
    x = x(:);
    xi = xi(:);
    N = numel(xi);
    L = zeros(N, 1);
    R = zeros(N, 1);
    for k = 1:N
        less = x <= xi(k);
        more = x >= xi(k);
        if any(less)
            L(k) = find(less, 1, 'last');
        end
        if any(more)
            R(k) = find(more, 1, 'first');
        end
    end
    u = L & R; % valid interpolatory points
    I = find(u);
    L = L(u);
    R = R(u);
    xi = xi(u);
    
    I = [I I];
    J = [L R];
    dx = x(R) - x(L);
    V = [(x(R)-xi), (xi-x(L))];
    V(dx>0) = spdiag( 1./dx(dx>0) ) * V(dx>0);
    V(~dx) = 1;
    S = sparse(I, J, V, N, numel(x));
end
