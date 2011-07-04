classdef Solution < handle

properties
    grid;
    vars;
end

methods
    function self = Solution(varargin)
        self.grid = struct(varargin{:}); 
        self.vars = struct();
        names = fieldnames(self.grid);
        n = 0;
        for k = 1:numel(names)
            g = self.grid.(names{k});
            n = n + g.numel;
        end
        self.grid.numel = n;
        offset = 0;
        for k = 1:numel(names)
            g = self.grid.(names{k});
            n = g.numel;            
            I = offset + (1:n);
            self.vars.(names{k}) = Selector(g, self, names{k}, I);
            offset = offset + n;
        end
    end
    function g = grad(self)
        g = 1;
    end
    function r = res(self)
        r = self.
    end
end

end

