function step = newton_solver1(grid)

    %% Differential operators 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Divergence, Gradient and Interpolation
    [D1, G1, I1] = operators(grid.center, 1);
    [D2, G2, I2] = operators(grid.center, 2);
    L = sparse_laplacian(grid.center);    
    
    function [rhs, H] = full_system(sol)
        % Compute gradients and interpolated values
        I1_C = I1 * sol.C(:); % C on X-staggered grid
        I2_C = I2 * sol.C(:); % C on Y-staggered grid
        G1_C = G1 * sol.C(:); % dC/dx on X-staggered grid
        G2_C = G2 * sol.C(:); % dC/dy on Y-staggered grid
        
        % Advection
        P1 = upwind(grid.Vx, sol.Vx, 1); % Choice of upwind cell's face
        P2 = upwind(grid.Vy, sol.Vy, 2); % Choice of upwind cell's face
        
        % Remove ghost points from V grid.
        R1 = P1 * select(grid.Vx.I | shift(grid.Vx.I, [1 0]) | shift(grid.Vx.I, [-1 0]));
        R2 = P2 * select(grid.Vy.I | shift(grid.Vy.I, [0 1]) | shift(grid.Vy.I, [0 -1]));
        % Select velocity from an upwind face
        R_Vx = R1 * sol.Vx(:);
        R_Vy = R2 * sol.Vy(:);
        % Select grad(C) from an upwind face
        PG_Cx = P1 * G1_C;
        PG_Cy = P2 * G2_C;
        % Advection upwind term
        V_gradC = R_Vx .* PG_Cx + R_Vy .* PG_Cy;

        % Right-Hand Side (operator applied to real- and ghost- values)        
        rhs = [D1 * G1_C + D2 * G2_C - sol.alpha * V_gradC];
           
        V_grad = spdiag(R_Vx) * P1 * G1 + spdiag(R_Vy) * P2 * G2;
         
        H = L - sol.alpha * V_grad; % C
    end

    %% Boundary conditions 
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Neumann on R=inf
    H_C0 = ident(grid.C);
    S_C = symmetries(grid.C);
    
    function [sol, H] = boundary_conditions(sol)

        sol.C(1, 2:end-1) = 1;
        sol.C(end, 2:end-1) = 0;
        sol.C(1:end, 1) = sol.C(1:end, 2);
        sol.C(1:end, end) = sol.C(1:end, end-1);
        % Dirichlet (1) on R=inf
        % Non-linear coupling on R=1
        % Symmetry on Theta=0 & pi.
        H = S_C * H_C0;

        J = grid.C.I(:);
        H = H * expand(J);
        
        sol.u = [sol.Phi(:)];
    end

    %% Applies step to system interior variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sol = apply(sol, du)
        sol.C(grid.C.I) = sol.C(grid.C.I) + du(:);
    end

    %% Newton method step
    %%%%%%%%%%%%%%%%%%%%%%%%
    function [sol, res, du] = newton_step(sol)
        [sol, Hb] = boundary_conditions(sol); % Apply boundary conditions
        [res, Hf] = full_system(sol); % Apply full system to get residual -> 0
        sol.res = res;
        % 0 = F(B(u)) + Hf * Hb * du
        % H * du = -F(B(u)) = -residual
        % du = - H \ residual
        H = Hf * Hb; % Interior's Hessian
        du = -(H \ res);
        sol = apply(sol, du);
    end
    step = @newton_step;
    
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
