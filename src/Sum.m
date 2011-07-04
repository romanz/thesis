classdef Sum < Binary

methods
    function self = Sum(op1, op2)
        self = self@Binary(op1, op2);
    end
    function r = res(self)
        r = self.op1.res() + self.op2.res();
    end
    function g = grad(self)
        g = self.op1.grad() + self.op2.grad();
    end
end

end