function [S, Kx, Ky] = stencil(s, Dx, Dy, Ix, Iy)

% We assume the dimensions are (x, y).
% ndgrid is the right function - not meshgrid.
% [Dx, Dy] = ndgrid(Dx, Dy);
% [Ix, Iy] = ndgrid(Ix, Iy);
% Dx = Dx(:);
% Dy = Dy(:);
% Ix = Ix(:);
% Iy = Iy(:);
Kx = repmat(Ix, [1 numel(Dx)]); 
Ky = repmat(Iy, [1 numel(Dy)]); 
Dx = repmat(Dx.', [numel(Ix) 1]);
Dy = repmat(Dy.', [numel(Iy) 1]);
S = laplacian(s, Dx, Dy, Kx, Ky);
Kx = Kx + Dx;
Ky = Ky + Dy;

function S = laplacian(s, Dx, Dy, varargin)
assert(all(size(Dx) == size(Dy)))
S = zeros(size(Dx));
J = (abs(Dx) <= 1) & (abs(Dy) <= 1);
S(J) = s(sub2ind(size(s), Dx(J)+2, Dy(J)+2));
