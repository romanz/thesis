function [M] = stokes_vanka_redblack(sz, A)
    sz = sz - 2;
    K = 1:prod(sz);
    K = K(:);
    [I, J] = ind2sub(sz, K);
    
    % Same indices for variables and equations.
    % CR: make this offset hack better.
    V = [[index(sz - [1 0], I-1, J), index(sz - [1 0], I, J)], ...
         [index(sz - [0 1], I, J-1), index(sz - [0 1], I, J)] + prod(sz - [1 0]), ...
         K + prod(sz - [1 0]) + prod(sz - [0 1])];
    
    K = logical(mod(I - J, 2)); 
    V1 = arr2cell(V( K, :)); % odd
    V2 = arr2cell(V(~K, :)); % even
    
    M{1} = vanka(A, V1, V1);
    M{2} = vanka(A, V2, V2);
end

function C = arr2cell(A)
    C = cell(size(A, 1), 1);
    for k = 1:numel(C)
        a = A(k, :);
        a = a(~isnan(a));
        C{k} = a(:);
    end
end

function K = index(sz, I, J)
    Q = (1 <= I) & (I <= sz(1)) & (1 <= J) & (J <= sz(2));
    K = nan(size(Q));
    K(Q) = sub2ind(sz, I(Q), J(Q));
end
