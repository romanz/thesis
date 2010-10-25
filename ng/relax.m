% Iterative relaxation
%   [U, RESIDUAL] = RELAX(M, A, F, U, ITERS)
function [u, res] = relax(M, A, f, u, iters)
    for iter = 1:iters
        for i = 1:numel(M)
            res = f - A*u;
            u = u + M{i}*res;
        end
    end
end
