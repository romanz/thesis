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
        else
            self.val = val(grid.R, grid.T);
        end
    end
    function r = res(self)
        r = self.val(:);
    end
    function G = grad(self)
        G = 0;
    end
end

end