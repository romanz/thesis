% Regular structured grid class
classdef Grid < handle
properties
    r, t;
    R, T;
    size;    
    numel;
end

methods
    function self = Grid(radius, theta)
        self.r = radius(:);
        self.t = theta(:);
        [self.R, self.T] = ndgrid(self.r, self.t);
        self.size = [numel(self.r), numel(self.t)]; %#ok<CPROP,PROP>
        self.numel = prod(self.size);
    end
    
    function res = eq(g1, g2)
        res = isequal(g1.r, g2.r) && isequal(g1.t, g2.t);
    end
    
    function g = crop(self, br, bt)
        gr = self.r(1+br:end-br);
        gt = self.t(1+bt:end-bt);
        g = Grid(gr, gt);
    end
    
    function g = boundary(self, b, dim)
        gr = self.r;
        gt = self.t;
        switch dim
            case 1, gr = [gr(1:b); gr(end+b+1:end)];
            case 2, gt = [gt(1:b); gt(end+b+1:end)];
        end
        g = Grid(gr, gt);
    end
end

end
