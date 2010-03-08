% Performs in-place QR using Modified Gram-Schmidt algorithm
%   [Q, R] = MGS(A), such that Q*R = A, Q'*Q = I
%   and R is upper-diagonal.
%   Note that we replace A by Q, thus saving memory.
%
% Roman Zeyde, Computer Science Department
% Technion -- Israel Institute of Technology
% romanz@cs.technion.ac.il

function [Q, R] = MGS(Q)

[N, M] = size(Q);
R = zeros(M);
for i = 1:M
    for j = 1:i-1
        r = Q(:, j)' * Q(:, i);
        Q(:, i) = Q(:, i) - r * Q(:, j);
        R(j, i) = r;
    end
    R(i, i) = norm(Q(:, i));
    Q(:, i) = Q(:, i) / R(i, i);
end
