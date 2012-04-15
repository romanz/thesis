classdef Minus < Binary

methods
    function self = Minus(op1, op2)
        self = self@Binary(op1, op2);
    end
    function r = res(self)
        r = self.op1.res() - self.op2.res();
    end
    function g = grad(self)
        g1 = self.op1.grad();
        g2 = self.op2.grad();
        g = g1 - g2;
    end
end

end
