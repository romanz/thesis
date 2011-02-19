function [D, G, I] = operators(grid, dim)
    dir = (1:2) == dim;
    assert(any(dir));

    J = grid.I;
    J = J | shift(J, dir) | shift(J, -dir);

    [G, I] = spdiff(J, dim);
    D = spdiff(true(grid.sz + dir - 2), dim);
    switch dim
        case 1, 
            a = grid.X;
            G = spdiag( 1./(G * a(:)) ) * G;
            D = spdiag( 1./(D*I*grid.X(:)) ) * D;
            D = D * spdiag( (I * grid.X(:)).^2 );
            D = spdiag( grid.X(grid.I).^(-2) ) * D;
        case 2, 
            a = grid.X .* grid.Y;
            G = spdiag( 1./(G * a(:)) ) * G;
            D = spdiag( 1./(D*I*grid.Y(:)) ) * D;
            D = D * spdiag( sin(I * grid.Y(:)) );
            D = spdiag( 1 ./ sin(grid.Y(grid.I)) ) * D;
            D = spdiag( grid.X(grid.I).^(-1) ) * D;
        otherwise
            error('invalid dimension');
    end
end

% Sparse matrix gradient representation for finite difference operator.
function [D, I] = spdiff(mask, dim)
    dir = (1:2 == dim); assert(any(dir));
    J1 = shift(mask, dir);
    J0 = shift(J1, -dir);
    J = [find(J0(:)) find(J1(:))]; % column indices
    N = size(J, 1); % # of equations
    M = numel(mask); % # of variables
    I = repmat(1:N, 1, 2); % row indices
    
    % Finite difference operator
    D = sparse(I, J, repmat([-1 1], N, 1), N, M); 
    I = sparse(I, J, repmat(1/2, N, 2), N, M); 
end
