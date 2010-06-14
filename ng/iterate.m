%% Solve problem using iterative process.
function varargout = iterate(prob, varargin)
    sz = zeros(nargin, 2);
    for i = 1:numel(varargin)
        sz(i, :) = size(varargin{i});
        varargin{i} = varargin{i}(:);
    end
    x = cat(1, varargin{:});
    for i = 1:prob.iters
        for j = 1:numel(prob.precond)
            residual = prob.rhs - prob.operator * x;
            dx = prob.precond{j} * residual;
            x = x + dx;
        end
    end
    k = 0;
    for i = 1:numel(varargin)
        n = prod(sz(i, :));
        varargout{i} = reshape(x(k + (1:n)), sz(i, :));
        k = k + n;
    end
end