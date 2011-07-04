classdef Deriv < Linear

methods
    function self = Deriv(grid, op, dim)
        self = self@Linear(grid, op);
        
        dir = (1:2) == dim;
        assert(any(dir));

        J = grid.I;
        J = J | shift(J, dir) | shift(J, -dir);

        J1 = shift(J, dir); % forward mask
        J0 = shift(J1, -dir); % backward mask
        J = [find(J0(:)) find(J1(:))]; % column indices
        
        N = size(J, 1); % # of equations
        M = numel(mask); % # of variables
        I = repmat(1:N, 1, 2); % row indices

        % Finite difference operator
        L = sparse(I, J, repmat([-1 1], N, 1), N, M); 
        
        switch dim % Choose coordinate
            case 1, t = grid.X;
            case 2, t = grid.Y;
        end
        
        L = spdiag( 1./(L * t(:)) ) * L;
        self.L = L;
    end
end

end