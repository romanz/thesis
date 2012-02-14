function L = scalar_laplacian(f)
    L = simplify(divergence(grad(f)));
end
