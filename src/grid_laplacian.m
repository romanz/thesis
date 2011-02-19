% Create Laplace operator L on grid G.
function L = grid_laplacian(G)
L = laplacian(G.I, G.X, G.Y);
