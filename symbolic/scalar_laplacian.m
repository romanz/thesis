function L = scalar_laplacian(f)
    L = simplify(divergence(gradient(f)));
end
