function [G, I] = sparse_gradient(grid, dim)
    dir = (1:2) == dim;
    assert(any(dir));

    J = grid.I;
    J = J | shift(J, dir) | shift(J, -dir);
    [G, I] = spdiff(J, dim);
    switch dim
        case 1, 
            a = grid.X;
            G = spdiag( 1./(G * a(:)) ) * G;
        case 2, 
            a = grid.X .* grid.Y;
            G = spdiag( 1./(G * a(:)) ) * G;
        otherwise
            error('invalid dimension');
    end
end
