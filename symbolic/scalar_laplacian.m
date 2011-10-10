function L = scalar_laplacian(f)
    L = divergence(gradient(f));
end
