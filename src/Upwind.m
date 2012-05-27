% Linear Sparse operator.
classdef Upwind < Operator
properties
    C, V;
    updater;
end

methods
    function self = Upwind(C, V)
        self = self@Operator(V.grid);
        self.C = C;
        self.V = V;
        
        dim = (C.grid.size - V.grid.size) > 0; % [1 0] for R, [0 1] for T.
        assert(sum(dim) == 1);
        
        I = 1:V.grid.numel;
        J = true(self.C.grid.size); % C grid
        J1 = find(shift(J, +dim)); % Positive cell
        J0 = find(shift(J, -dim)); % Negative cell
        I = [I, I];
        J = [J0, J1];
        spinit = @(I, J, sz) (@(val) sparse(I, J, val, sz(1), sz(2)));
        self.updater = spinit(I, J, [self.V.grid.numel, self.C.grid.numel]);
    end
    function L = op(self)
        v = self.V.res();
        L = self.updater([v>0, v<0]); 
        % Take upstream cells:
        % * Negative cells for positive velocity.
        % * Positive cells for negative velocity.
    end
    function r = res(self)
        r = self.op() * self.C.res();
    end
    function G = grad(self)
        G = dot_prod(self.op(), self.C.grad());
    end
end

end