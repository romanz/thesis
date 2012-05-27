function [D] = sparse_divergence(grid, dim, central)
    dir = (1:2) == dim;
    assert(any(dir));
    J = grid.I;
    J = J | shift(J, dir) | shift(J, -dir);
    D = spdiff(J, dim);
    cX = central.X(central.I);
    cY = central.Y(central.I);
    switch dim
        case 1,
            D = spdiag( 1./(D*grid.X(:)) ) * D; % d(Fr)/dr
            D = D * spdiag( grid.X(:).^2 ); % d(r^2 Fr)/dr
            D = spdiag( cX(:).^(-2) ) * D; 
            % 1/r^2 d(r^2 Fr)/dr
        case 2, 
            D = spdiag( 1./(D*grid.Y(:)) ) * D; % d(Ft)/dt
            D = D * spdiag( sin( grid.Y(:)) ); % d(sint Ft)/dt
            D = spdiag( 1 ./ (cX(:) .* sin(cY(:))) ) * D; 
            % 1/(r sint) d(sint Ft)/dt
        otherwise
            error('invalid dimension');
    end
end
