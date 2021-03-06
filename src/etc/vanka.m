function [M] = vanka(A, V, E)
% Vanka-type smoother construction.

% V is a cell array, whose k-th entry contains k-th subdomain indices
nonzeros = 0;
for k = 1:numel(V)
    nonzeros = nonzeros + numel(V{k}).^2;
end
% Preallocate memory for non-zeroes
inverses = zeros(nonzeros, 1);
indices = zeros(size(inverses));

sz = size(A); % Pre-compute A size
offset = 0;
for k = 1:numel(V)
     % Process current index set:
    J = col(V{k});
    I = col(E{k});
    N = numel(J); 
    % Invert current submatrix
    inverses( offset + (1:N^2) ) = inv( A(I, J) );
    % Save its indices at A
    indices( offset + (1:N^2) ) = index_ndgrid(sz, J, I);
    offset = offset + N^2;
end
% Construct sparse matrix M for pre-conditioning
M = mksparse(sz, indices, inverses);

% r = f - Ax
% x' = x + M r 
% x' = (I - MA)x + Mf
% x' = Tx + d

function K = index_ndgrid(sz, I, J)
[I, J] = ndgrid(I, J);
K = sub2ind(sz, I, J);

function S = mksparse(sz, indices, values)
[I, J] = ind2sub(sz, indices);
S = sparse(I, J, values, sz(1), sz(2));
