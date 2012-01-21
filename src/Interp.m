classdef Interp < Linear

methods
    function self = Interp(grid, op)
        L = interpolate_grids(grid, op.grid);
        self = self@Linear(grid, op, L);
    end
end

end

function S = interpolate_grids(dst, src)
    n = fieldnames(src);
    S1 = interpolate_matrix(src.(n{1}), dst.(n{1}));
    S1 = kron(speye(numel(src.(n{2}))), S1); % Interpolate on 1st dimension
    S2 = interpolate_matrix(src.(n{2}), dst.(n{2}));
    S2 = kron(S2, speye(numel(dst.(n{1})))); % Interpolate on 2nd dimension
    S = S2 * S1;    
end

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
    V(dx>0, :) = spdiag( 1./dx(dx>0) ) * V(dx>0, :);
    V(~dx, :) = 1/2;
    S = sparse(I, J, V, N, numel(x));
end
