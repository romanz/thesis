classdef Interp < Linear

methods
    function self = Interp(grid, op)
        self = self@Linear(grid, op);
        self.L = interpolate_grids(op.grid, grid);
    end
end

end

function S = interpolate_grids(src, dst)
    n = fieldnames(src);
    S1 = interpolate_matrix(src.(n{1}), dst.(n{1}));
    S1 = kron(speye(numel(src.(n{2}))), S1); % Interpolate on 1st dimension
    S2 = interpolate_matrix(src.(n{2}), dst.(n{2}));
    S2 = kron(S2, speye(numel(dst.(n{1})))); % Interpolate on 2nd dimension
    S = S2 * S1;    
end
