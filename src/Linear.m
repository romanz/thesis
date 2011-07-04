classdef Linear < Operator
properties
    L;
    op;
end

methods
    function self = Linear(grid, op)
        self = self@Operator(grid);
        self.op = op;
    end
    function r = res(self)
        r = self.L * self.op.res();
    end
    function G = grad(self)
        G = self.L * self.op.grad();
    end
end

end