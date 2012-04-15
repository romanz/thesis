classdef Selector < Interp

methods
    function self = Selector(grid, op)
        self = self@Interp(grid, op);
        assert(all(nonzeros(self.L) == 1));
    end
end

end

