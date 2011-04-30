function step = newton_solver(grid)

    function [sol, rhs, du] = newton_step(sol)
        [sol, Hb] = boundary_conditions(sol); % Apply boundary conditions
        [rhs, Hf] = full_system(sol); % Apply full system to get RHS -> 0
        % 0 = F(B(u)) + Hf * Hb * du
        % H * du = -F(B(u)) = -rhs
        % du = - H \ rhs
        H = Hf * Hb; % Interior's Hessian
        du = -(H \ rhs);
        sol = apply(sol, du);
    end
    step = @newton_step;
    
    % Divergence, Gradient and Interpolation
    [D1, G1, I1] = operators(grid.center, 1);
    [D2, G2, I2] = operators(grid.center, 2);
    L = sparse_laplacian(grid.center);
    
    function [rhs, H] = full_system(sol)
        % Compute gradients and interpolated values
        I1_C = I1 * sol.C(:); % C on X-staggered grid
        I2_C = I2 * sol.C(:); % C on Y-staggered grid
        G1_Phi = G1 * sol.Phi(:); % dPhi/dx on X-staggered grid
        G2_Phi = G2 * sol.Phi(:); % dPhi/dy on Y-staggered grid
        G1_C = G1 * sol.C(:); % dC/dx on X-staggered grid
        G2_C = G2 * sol.C(:); % dC/dy on Y-staggered grid
        
        % Right-Hand Side (operator applied to real- and ghost- values)        
        rhs = [D1 * (I1_C .* G1_Phi) + D2 * (I2_C .* G2_Phi); ...
               D1 * G1_C + D2 * G2_C];
           
        H11 = D1 * spdiag(I1_C) * G1 + D2 * spdiag(I2_C) * G2; % Phi
        H12 = D1 * spdiag(G1_Phi) * I1 + D2 * spdiag(G2_Phi) * I2; % C
        
        H21 = sparse(size(L, 1), grid.Phi.numel); % Phi
        H22 = L; % C
        
        H = [H11 H12; H21 H22];
    end

    % Neumann on R=inf
    H_Phi0 = ident(grid.Phi) + coupling(grid.Phi, [1 0]);
    H_C0 = ident(grid.C);
    S_Phi = symmetries(grid.Phi);
    S_C = symmetries(grid.C);

    function [sol, H] = boundary_conditions(sol)

        Er = -sol.beta * cos(grid.Phi.y(2:end-1)).';
        dr = grid.Phi.x(end) - grid.Phi.x(end-1);

        Phi = sol.Phi;
        Phi0 = -log( sol.C(2, 2:end-1) ); % Phi = -log(C)
        dPhi0_dC = -1 ./ sol.C(2, 2:end-1); % dPhi/dC = -1/C
        Phi(1, 2:end-1) = Phi0;
        Phi(end, 2:end-1) = Phi(end-1, 2:end-1) + Er * dr;
        Phi(1:end, 1) = Phi(1:end, 2);
        Phi(1:end, end) = Phi(1:end, end-1);
        % Non-linear coupling on R=1
        % Symmetry on Theta=0 & pi.
        H_Phi = S_Phi * ([H_Phi0, coupling(grid.C, [-1 0], dPhi0_dC)]);

        C = sol.C;
        C0 = exp( -sol.Phi(2, 2:end-1) ); % C = exp(-Phi)
        dC0_dPhi = -C0; % dC/dPhi = -exp(-Phi) = -C
        C(1, 2:end-1) = C0;
        C(end, 2:end-1) = 1;
        C(1:end, 1) = C(1:end, 2);
        C(1:end, end) = C(1:end, end-1);
        % Dirichlet (1) on R=inf
        % Non-linear coupling on R=1
        % Symmetry on Theta=0 & pi.
        H_C = S_C * ([coupling(grid.Phi, [-1 0], dC0_dPhi), H_C0]);

        H = [H_Phi; H_C];
        
        H = H * expand([grid.Phi.I(:); grid.C.I(:)]);

        sol.Phi = Phi;
        sol.C = C;
    end

    function sol = apply(sol, du)
        [dPhi, dC] = split(du, ...
            grid.Phi.sz-2, grid.C.sz-2);
        sol.Phi(grid.Phi.I) = sol.Phi(grid.Phi.I) + dPhi(:);
        sol.C(grid.C.I) = sol.C(grid.C.I) + dC(:);
    end
end

function x = remove_mean(x)
    x = x - mean(x(:));
end

function S = ghost(I, dir, H)
    if nargin < 3
        H = 1;
    else
        H = spdiag(H);
    end
    J = ~I & shift(I, dir); % Ghost points
    S = expand(J) * H * select(shift(J, -dir));
end

function S = symmetries(grid)
    J = true(grid.sz);
    J(:, [1 end]) = false;
    S = expand(J) * select(J) + ghost(J, [0 -1]) +  ghost(J, [0 1]);
end

function S = coupling(grid, varargin)
    S = ghost(grid.I, varargin{:});
end

function S = ident(grid)
    K = find(grid.I(:));
    S = sparse(K, K, 1, grid.numel, grid.numel);
end

% Splits x into seperate variables according to specified sizes.
%   [x1, x2, x3] = split(x, sz1, sz2, sz3);
%   where sz{i} = size(x{i}) and x = [x1(:); x2(:); x3(:)];
function [varargout] = split(x, varargin)
    N = numel(varargin);
    varargout = cell(N);
    offset = 0;
    for k = 1:N
        sz = varargin{k};
        m = prod(sz);
        varargout{k} = reshape(x(offset + (1:m)), sz);
        offset = offset + m;
    end
    assert(offset == numel(x), 'Splitting mismatch.');
end
