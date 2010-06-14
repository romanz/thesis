function [lhs, rhs] = neumann(grid, dir)
    J = ~grid.I & shift(grid.I, dir); % boundary logical mask
    I = shift(J, -dir); % interior neighbourhood logical mask
    Ji = find(J); % boundary indices
    Ii = find(I); % interior indices
    h = [grid.X(J) - grid.X(I), grid.Y(J) - grid.Y(I)];
    h = h * dir(:);
    N = prod(grid.sz);
    K = (1:N)';
    lhs = sparse([K; Ji], [K; Ii], [double(~J(:)); ones(size(Ji))], N, N); % A = A*P
    rhs = @func;
    function q = func(vals)
        if numel(vals) == 1
            vals = repmat(vals, [numel(nnz(J), 1) 1]);
        end
        q = sparse(Ji, ones(size(Ji)), -h(:) .* vals(:), N, 1); % b = A*q
    end
end
