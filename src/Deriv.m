classdef Deriv < Linear

methods
    function self = Deriv(grid, op, dim)
        self = self@Linear(grid, op);
        
        dir = (1:2) == dim;
        assert(any(dir));

        J = true(op.grid.sz);
        J1 = shift(J, dir); % forward mask
        J0 = shift(J1, -dir); % backward mask
        J = [find(J0(:)) find(J1(:))]; % column indices
        
        N = size(J, 1); % # of equations
        M = op.grid.numel; % # of variables
        I = repmat(1:N, 1, 2); % row indices

        % Finite difference operator
        L = sparse(I, J, repmat([-1 1], N, 1), N, M); 
        
        switch dim % Choose coordinate
            case 1, a = op.grid.R;
            case 2, a = op.grid.T;
        end
        
        L = spdiag( 1./(L * a(:)) ) * L;
        self.L = L;
    end
end

end