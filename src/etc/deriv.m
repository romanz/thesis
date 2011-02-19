% D = DERIV(G, DIR)
%   Derivation operator on a grid on coordinates, in specified direction.
%   Example:
%       [X, Y] = ndgrid(x, y);
%       Dx = deriv(X, [1 0]);
%       Dy = deriv(Y, [0 1]);
%
function D = deriv(G, dir)
I = true(size(G));
I1 = shift(I, -dir);
I2 = shift(I, +dir);
D = 1 ./ (G(I2) - G(I1));
N = prod(size(G)-dir);
D = sparse(repmat((1:N)', [1 2]), [find(I1), find(I2)], [-D, D]);
