classdef Deriv < Linear

methods
    function self = Deriv(grid, op, dim)
        dir = (1:2) == dim;
        assert(any(dir));

        J = true(op.grid.size);
        J1 = shift(J, dir); % forward mask
        J0 = shift(J1, -dir); % backward mask
        J = [find(J0(:)) find(J1(:))]; % column indices
        
        N = size(J, 1); % # of equations
        M = op.grid.numel; % # of variables
        I = repmat(1:N, 1, 2); % row indices

        % Finite difference operator
        D = sparse(I, J, repmat([-1 1], N, 1), N, M); 
        
        switch dim % Choose coordinate
            case 1, x = op.grid.R;
            case 2, x = op.grid.T;
        end
        
        % L = D/Dx
        Dx = D * x(:);
        assert(all(Dx));
        
        L = spdiag( 1./Dx ) * D;
        self = self@Linear(grid, op, L);        
    end
end

end
