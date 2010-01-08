% Create Jacobi iteration matrix B = inv(diag(A))
% [B] = JACOBI(A)
function [B] = jacobi(A)
% A = D + LU
% f = Av = Dv + LUv
% v' = D^{-1}(f - LUv)
% r = f - Av = f - Dv - LUv
% v' = D^{-1}(r + Dv) = v + D^{-1}r

if nargin < 2
    w = 1;
end

N = length(A);
I = sub2ind(size(A), 1:N, 1:N);
I = I(:);

D = A(I);
B = sparse(1:N, 1:N, 1./D);
