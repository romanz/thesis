classdef Join < handle
properties
    ops;
end

methods
    function self = Join(varargin)
        self.ops = varargin;
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
