function [R, g] = jacobi(A, f, w)
if nargin < 3
    w = 1;
end

N = numel(f);
I = sub2ind(size(A), 1:N, 1:N);
I = I(:);

D = A(I);
A(I) = 0;
R = -A .* repmat(1./D, [1 N]);

R = (1 - w) * eye(N) + w * R;
g = f ./ D;
