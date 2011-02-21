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
