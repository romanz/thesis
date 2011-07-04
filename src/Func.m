classdef Func < Operator
properties
    op; % Input operator
    func; % y = f(x)
    deriv; % dy/dx = f'(x)
end

methods
    function self = Func(op, func)
        grid = op.grid; % Reuse operator's grid
        self = self@Operator(grid);
        self.op = op;
        if ~isempty(func)
            if isa(func, 'char')
                func = sym(func);
            end % using symoblic toolbox
            self.func = matlabFunction(func);
            self.deriv = matlabFunction(diff(func));
        else % just copy values to grid
            self.func = @(x) x;
            self.deriv = @(x) ones(size(x));
        end
    end
    function r = res(self)
        r = self.func( self.op.res() );
    end
    function G = grad(self)
        r = self.op.res();
        G = spdiag( self.deriv(r) ) * self.op.grad();
    end
end

end