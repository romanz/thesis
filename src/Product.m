classdef Product < Binary

methods
    function self = Product(op1, op2)
        self = self@Binary(op1, op2);
    end
    function r = res(self)
        r = self.op1.res() .* self.op2.res();
    end
    function g = grad(self)
        g = spdiag(self.op2.res()) * self.op1.grad() + ...
            spdiag(self.op1.res()) * self.op2.grad();
    end
end

end