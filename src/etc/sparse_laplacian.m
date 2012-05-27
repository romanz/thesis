function L = sparse_laplacian(grid)
    [D1, G1] = operators(grid, 1);
    [D2, G2] = operators(grid, 2);
    L = D1 * G1 + D2 * G2;
end
