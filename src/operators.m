function [D, G, I] = operators(grid, dim)
    dir = (1:2) == dim;
    assert(any(dir));

    [G, I] = sparse_gradient(grid, dim);

    switch dim
        case 1, hx = [1;1]/2; hy = [0;1;0];
        case 2, hy = [1;1]/2; hx = [0;1;0];
    end
    x = convn(grid.x, hx, 'valid');
    y = convn(grid.y, hy, 'valid');
    D = sparse_divergence(init_grid(x, y, false), dim, grid);    
end
