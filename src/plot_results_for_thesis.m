function plot_results_for_thesis(filename)

% beta = 3;
% N = 512;
% Rmax = 100;
% filename = sprintf('~/MATs/sol_beta=%.3e_[%dx%d]_Rmax=%.1f_Du=1.00_zeta=10.00_alpha=0.50', ...
%             beta, N+1, N+1, Rmax);
% load([filename '.mat'])
load(filename)
[p, n] = fileparts(filename);
filename = fullfile(p, n);
disp(filename)
if 1
    clf;
    hold on;
    op = sol.C;
    Z = regrid(op);
    g = op.grid;
    h = [1 1]/2;
    R = convn(g.R, h, 'valid');
    T = convn(g.T, h, 'valid');
    z = convn(Z, h, 'valid');
    x = R .* cos(T);
    y = R .* sin(T);
    c = 1 + linspace(-1, 1, 101) * sol.beta;
    contourf(x, y, z, c, '-');
    contourf(x, -y, z, c, '-');
    t = linspace(0, 2*pi, 100);
    fill(cos(t), sin(t), [1 1 1]*0.5)
    axis([-1 1 -1 1]*2.5)
    axis equal
    colorbar
    title(sprintf('Salt concentration (C): \\beta = %.2f', sol.beta))
    print('-depsc2', [filename '_C.eps'])
end
if 0
    for d = [1 -1]
        clf;
        hold on;
        s = streamfunc(sol);
        I = find(s.grid.Psi.r < 10);
        R = s.grid.Psi.R(I, :);
        T = s.grid.Psi.T(I, :);
        z = s.Psi(I, :);
        x = R .* cos(T);
        y = R .* sin(T);
        c = linspace(-1, 1, 101)*sol.beta*10;
        hold on;
        contour(d*x, y, z, c)
        contour(d*x, -y, z, c)
        t = linspace(0, 2*pi, 100);
        fill(cos(t), sin(t), [1 1 1]*0.5)
        axis([-1 1 -1 1]*3)
        axis equal
        title(sprintf('Streamlines (\\Psi): \\beta = %.2f', sol.beta))
        print('-depsc2', [filename sprintf('_[%d]_Psi.eps', d)])
    end
end

if 0
    clf;
    hold on;
    op = sol.Phi;
    g = op.grid;
    I = g.r < 6;
    R = g.R(I, :);
    T = g.T(I, :);
    X = R .* cos(T);
    Y = R .* sin(T);
    Z = regrid(op);
    Z = Z(I, :);
    Er = diff(Z(1:2, :), 1, 1) ./ diff(g.r(1:2));
    Er = (Er(1:end-1) + Er(2:end)) / 2;
    Er = -Er(:);
    Et = diff(mean(Z(1:2, :)), 1, 2) ./ diff(g.t(1:2));
    Et = -Et(:);
    t = sol.grid.t;
    Ex = Er .* cos(t) - Et .* sin(t);
    Ey = Er .* sin(t) + Et .* cos(t);
    x = convn(X, [1 1; 1 1]/4, 'valid');
    y = convn(Y, [1 1; 1 1]/4, 'valid');
    z = convn(Z, [1 1; 1 1]/4, 'valid');
    z = z - z(1, end);
    c = (-1:0.003:1)*sol.beta*10;
    contour(x, y, z, c)
    contour(x, -y, z, c)
    t = linspace(0, 2*pi, 1000);
    fill(cos(t), sin(t), [1 1 1]*0.5)
    x = cos(sol.grid.t); y = sin(sol.grid.t);
    scale = 0.1 / sol.beta;
    quiver([x; x], [y; -y], scale*[Ex; Ex], scale*[Ey; -Ey], 0, 'k')
    axis([-1 1 -1 1]*2)
    axis equal
    title(sprintf('Electric potential (\\phi) and field at the surface (E): \\beta = %.2f', sol.beta))
    print('-depsc2', [filename '_PhiE.eps'])
end

if 1
    sol = rmfield(sol, {'var', 'C', 'Phi', 'Vr', 'Vt', 'P', 'bnd', 'eqn', 'grid'});
    save([filename '_lite.mat'], 'sol')
end