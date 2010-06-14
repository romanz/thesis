function [lhs, rhs] = dirichlet(grid, dir)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    Ji = find(J); % boundary indices
    
    N = prod(grid.sz);
    lhs = sparse(1:N, 1:N, double(~J)); % A = A*P
    
    rhs = @func;
    function q = func(vals)
        if numel(vals) == 1
            vals = repmat(vals, [nnz(J) 1]);
        end
        q = sparse(Ji, ones(size(Ji)), -vals(:), N, 1); % b = A*q
    end
end
