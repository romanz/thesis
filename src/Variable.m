classdef Variable < handle
properties
    value;
    grid; % pseudo-grid (used only for .numel)
end

methods
    function self = Variable(varargin)
        self.value = col(varargin{:});
        self.grid.numel = numel(self.value);
    end
    
    function r = res(self)
        r = self.value(:);
    end
    
    function G = grad(self) 
        G = speye(numel(self.value));
    end
    
    function update(self, delta)
        assert(numel(self.value) == numel(delta))
        self.value = self.value + delta;
    end
end

end
