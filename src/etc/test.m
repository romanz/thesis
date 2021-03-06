G = grid([1 2 4 8], [0 1 2 3 4]);
z = round(10*randn(2, 3))
I = speye(G.numel);
[P1, Q1] = average_dirichlet(G, [-1 0]);
[P2, Q2] = neumann(G, [+1 0]);
R = I + neumann(G, [0 -1], true) + neumann(G, [0 +1], true);
s1 = [10;11;12];
s2 = [-1;-2;-3];
reshape(R*((I + P1 + P2)*expand(G.I)*z(:) + Q1*s1 + Q2*s2), G.sz)

% function test
%% Test spherical laplacian
% n = 64; 
% [X, Y] = ndgrid(logspace(log10(1), log10(4), n), linspace(0, pi, n)); I = interior([n n]);
% A = laplacian(I, X, Y);
% u = X .* cos(Y); u = u(:);
% v = A * u;
% e = norm(v, inf);
% assert(e < 1.25e-3)

%% Test spherical divergence
% n = 2^6;
% x = logspace(0, 1, n);
% y = linspace(0, pi, n);
% [center, interior, xstag, ystag] = grids(x, y);
% 
% gridV = xstag; dir = [1 0]; 
% K = gridV.I | shift(gridV.I, -dir);
% div_x = divergence(K, gridV.X, gridV.Y, dir);
% Vx = -cos(gridV.Y);
% 
% gridV = ystag; dir = [0 1]; 
% K = gridV.I | shift(gridV.I, -dir);
% div_y = divergence(K, gridV.X, gridV.Y, dir);
% Vy = sin(gridV.Y);
% 
% d = div_x * Vx(:) + div_y * Vy(:);
% e = norm(d(:), inf);
% mesh(reshape(d, [n n]-1))
% e
% assert(e < 1e-3)

%% Coupled simulation main script.
% 
% % Create grids
%     nx = 10;
%     ny = 10;
%     x = logspace(log10(1), log10(10), nx);
%     y = linspace(0, pi, ny);
%     [center, interior, xstag, ystag] = grids(x, y);
% 
% % Initialize variables
%     Vx = 0.*randn(xstag.sz);
%     Vy = 1+0.*randn(ystag.sz);
%     P = 0.*randn(interior.sz);
%     
%     problem = @(iters) struct('iters', iters);
% % Problems:
%     stokes = problem(0*10e3);
% 
% % Stokes problem
%     vel = 1;
% 
%     gridP = interior;
%     gridVx = xstag;
%     gridVy = ystag;
% 
%     % Spherical
%     [Gx_p, L_vx, Gx_vx] = stokes1(1, gridVx, gridP); % X dimension
%     [Gy_p, L_vy, Gy_vy] = stokes1(2, gridVy, gridP); % Y dimension
% 
%     [lhs1x, rhs1x] = dirichlet(gridVx, [-1 0]);
%     [lhs2x, rhs2x] = dirichlet(gridVx, [+1 0]);
%     [lhs3x] = neumann(gridVx, [0 -1]);
%     [lhs4x] = neumann(gridVx, [0 +1]);
%     lhsVx = lhs1x * lhs2x * lhs3x * lhs4x * restrict(gridVx.I);
%     Vr = -vel*cos(gridVx.y(2:end-1));
%     rhsVx = rhs1x(Vr) + rhs2x(Vr); 
% 
%     % average_
%     [lhs1y, rhs1y] = dirichlet(gridVy, [-1 0]);
%     [lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
%     [lhs3y] = dirichlet(gridVy, [0 -1]);
%     [lhs4y] = dirichlet(gridVy, [0 +1]);
%     lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
% 
%     L_vx1 = L_vx * lhsVx;
%     L_vy1 = L_vy * lhsVy;
%     Gx_vx1 = Gx_vx * lhsVx;
%     Gy_vy1 = Gy_vy * lhsVy;
%     
%     stokes.operator = [... 
%         L_vx1, sparse(size(L_vx1, 1), size(L_vy1, 2)), -Gx_p; ...
%         sparse(size(L_vy1, 1), size(L_vx1, 2)), L_vy1, -Gy_p; ...
%         Gx_vx1, Gy_vy1, sparse(prod(gridP.sz), prod(gridP.sz))];
%     stokes.precond = stokes_vanka_redblack(gridP.sz, stokes.operator);
%     
%     function update_stokes
%         % P (no boundary condition: 1 DOF, arbirary mean)
%         Vt = vel*sin(gridVy.y(2:end-1));
%         rhsVy = rhs1y(Vt) + rhs2y(Vt);
%         % 
%         
%         stokes.rhs = [...
%             L_vx * rhsVx(:); ...
%             L_vy * rhsVy(:); ...
%             Gx_vx * rhsVx(:) + Gy_vy * rhsVy(:)];
%         
%         stokes.expandVx = @(Vx) ...
%             reshape(lhsVx * Vx(:) - rhsVx, gridVx.sz);
%         stokes.expandVy = @(Vy) ...
%             reshape(lhsVy * Vy(:) - rhsVy, gridVy.sz);
%     end
% 
%     tic;
% 
%     % Iterate on problems
%         update_stokes;
%         [Vx1, Vy1, P1] = iterate(...
%             stokes, Vx(2:end-1, 2:end-1), Vy(2:end-1, 2:end-1), P);
%         [Vx1] = stokes.expandVx(Vx1);
%         [Vy1] = stokes.expandVy(Vy1);
%         P1 = P1 - mean(P1(:));
%         Vx = Vx1;
%         Vy = Vy1;
%         P = P1;
% 
%     toc;
%     figure(1); clf; show(gridVx, Vx(:, :), 'V_R')
%     figure(2); clf; show(gridVy, Vy(:, :), 'V_\theta')
%     figure(3); clf; show(gridP, P, 'P')
%     save results
% end

%% Vector laplacian test in Spherical coordinates
% clear
% fprintf('-------------------------------------------------------\n')
% 
% syms r t a
% Vr = cos(t)*(1 - 3*a/(2*r) + a^3/(2*r^3));
% Vt = -sin(t)*(1 - 3*a/(4*r) - a^3/(4*r^3));
% 
% subs(Vr, r, a)
% subs(Vr, r, inf)
% 
% subs(Vt, r, a)
% subs(Vt, r, inf)
% 
% Lr = diff(r^2 * diff(Vr, r), r) / r^2 + ...
%      diff(sin(t) * diff(Vr, t), t) / (r^2 * sin(t)) - ...
%      2 * (Vr + diff(sin(t) * Vt, t) / sin(t)) / r^2;
% Lr = simple(Lr)
% 
% Lt = diff(r^2 * diff(Vt, r), r) / r^2 + ...
%      diff(sin(t) * diff(Vt, t), t) / (r^2 * sin(t)) ...
%      - Vt / (r * sin(t))^2 + 2*diff(Vr, t) / r^2;
% Lt = simple(Lt)
% 
% P = -3*a/(2*r^2) * cos(t);
% 
% Fr = diff(P, r);
% Fr = simple(Fr)
% 
% Ft = diff(P, t) / r;
% Ft = simple(Ft)
% 
% (Lr == Fr && Lt == Ft)
% 
% div = diff(r^2 * Vr, r) / r^2 + diff(sin(t) * Vt, t) / (r * sin(t));
% div = simple(div)
% 
% Srr = 2*diff(Vr, r) - P; Srr = subs(Srr, r, a)
% Srt = diff(Vr, t) / r + r * diff(Vt/r, r); Srt = subs(Srt, r, a)
% Fz = Srr*cos(t) - Srt*sin(t);
% Fz = simple(Fz * 4*pi*a^2)
