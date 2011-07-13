classdef Linear < Operator
properties
    L;
    op;
end

methods
    function self = Linear(grid, op, L)
        self = self@Operator(grid);
        self.op = op;
        if isempty(L)
            L = speye(grid.numel); % Default copy operator
        end
        assert(size(L, 1) == grid.numel);
%         assert(size(L, 2) == op.grid.numel);
        self.L = L;
    end
    function r = res(self)
        r = self.L * self.op.res();
    end
    function G = grad(self)
        G = dot_prod(self.L, self.op.grad());
    end
end

end