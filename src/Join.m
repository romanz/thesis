classdef Join < handle
properties
    ops;
    grid; % pseudo-grid (used only for .numel)
end

methods
    function self = Join(varargin)
        self.ops = varargin;
        n = 0;
        for k = 1:numel(self.ops)
            op = self.ops{k};
            n = n + op.grid.numel;
        end
        self.grid.numel = n;
    end
    
    function r = res(self)
        r = cell(numel(self.ops), 1);
        for k = 1:numel(self.ops)
            r{k} = self.ops{k}.res();
        end
        r = cat(1, r{:});
    end
    
    function G = grad(self) 
        G = cell(numel(self.ops), 1);
        for k = 1:numel(self.ops)
            G{k} = self.ops{k}.grad();
        end
        G = cat(1, G{:});
    end
end
    
end
