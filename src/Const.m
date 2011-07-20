classdef Const < Operator
properties
    val;
end

methods
    function self = Const(grid, val)
        self = self@Operator(grid);
        if isa(val, 'float')
            if numel(val) == 1
                val = repmat(val, grid.sz);
            end
            self.val = val;
        else % function handle evaluated point-wise
            if isa(val, 'char')
                val = matlabFunction(val);
            end
            self.val = arrayfun(val, grid.R, grid.T);
        end
    end
    function r = res(self)
        r = self.val(:);
    end
    function G = grad(self)
        G = sparse(0); % [0]
    end
end

end