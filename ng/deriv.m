function D = deriv(G, dir)
I = true(size(G));
I1 = shift(I, -dir);
I2 = shift(I, +dir);
D = 1 ./ (G(I2) - G(I1));
N = prod(size(G)-dir);
D = sparse((1:N)'*[1 1], [find(I1), find(I2)], [-D, D]);
