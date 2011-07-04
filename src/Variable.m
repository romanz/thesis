classdef Variable < handle
properties
    value;
end

methods
    function self = Variable(varargin)
        self.value = col(varargin{:});
    end
    
    function r = res(self)
        r = self.value(:);
    end
    
    function G = grad(self) %#ok<MANU>
        G = 1; % speye(numel(self.value))
    end
end

end
