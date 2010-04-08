function [M, T, d] = vanka(A, f, S)
% Vanka-type smoother construction.

% S is a cell array, whose k-th entry contains k-th subdomain indices
nonzeros = 0;
for k = 1:numel(S)
    nonzeros = nonzeros + numel(S{k}).^2;
end
% Preallocate memory for non-zeroes
inverses = zeros(nonzeros, 1);
indices = zeros(size(inverses));

sz = size(A); % Pre-compute A size
offset = 0;
for k = 1:numel(S)
     % Process current index set:
    I = col(S{k});
    N = numel(I); 
    % Invert current minor
    inverses( offset + (1:N^2) ) = inv( A(I, I) );
    % Save its indices at A
    indices( offset + (1:N^2) ) = index_ndgrid(sz, I, I);
    offset = offset + N^2;
end
% Construct sparse matrix M for pre-conditioning
M = mksparse(sz, indices, inverses);

% r = f - Ax
% x' = x + M r 
% x' = (I - MA)x + Mf
% x' = Tx + d
T = speye(sz) - M * A;
d = M * f;

function K = index_ndgrid(sz, I, J)
[I, J] = ndgrid(I, J);
K = sub2ind(sz, I, J);

function S = mksparse(sz, indices, values)
[I, J] = ind2sub(sz, indices);
S = sparse(I, J, values, sz(1), sz(2));
