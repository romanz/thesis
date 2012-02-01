function show(name)
    load(name);
    figure(1); clf;
%     set(gcf, 'PaperPositionMode', 'auto');
%     set(gcf, 'Position',[100 100 900 600]);
    v = attr(solutions, 'Vinf');
    betas = attr(solutions, 'beta');

    g = sol.gamma;
    B3 = (31./(320*(g.^(1./2) + 1)) - 9./(320*(g.^(1./2) + 1).^2) + 1./1680);
    B1 = 2*log(1./(2*g.^(1./2)) + 1./2);
    b = logspace(log10(min(betas)), log10(max(betas)), 1000);
    w = b.^3*(B3 - B1*11./320) + b*B1;
    w3 = b.^3*(B3 - B1*11./320);
    w1 = b*B1;
    
    clf; 
    h = [];
    h = [h; graph(betas, v, '.')]; hold on; 
    h = [h; graph(b, w1+w3, '-')]; 
    graph(b, w1, ':'); 
    graph(b, w3, ':'); 
    title(sprintf('grid = [%d x %d], Rmax = %.1e, gamma = %.4f, alpha = %.4f', ...
        numel(sol.grid.radius), numel(sol.grid.theta), max(sol.grid.radius), ...
        sol.gamma, sol.alpha))
    xlabel('beta (electric field)')
    ylabel('U (steady-state velocity)')
    hold off;
    legend(h, {'Numerical (+)', 'Numerical (-)', ...
        'Theoretical (+)', 'Theoretical (-)'}, 'Location', 'Best');
    grid on; grid minor;
    print('-depsc2', [name(1:end-4) '.eps'])

    %{
    clf; hold on   
    for k = 1:numel(betas)
        s = solutions{k};
        c = mean(s.C(1:2, :));
        c = (((c - 1) ./ s.beta) + 1);
        phi = mean(s.Phi(1:2, :)) ./ s.beta;
        v = mean(s.Vy(1:2, :) ./ s.beta);
        plot(180 * s.grid.C.y ./ pi, [c; phi], 180 * s.grid.Vy.y ./ pi, v);
        xlim([0 180])
    end
    hold off;
    for k = 1:numel(betas)
        s = solutions{k};
        clf; hold on        
        % apply(@contour, s, 'Psi', 10, +1, s.Vinf * linspace(-1, 1, 50));
        b = s.beta;
        r = s.grid.Psi.X;
        t = s.grid.Psi.Y;
        g = s.gamma;
        psi = (b*sin(t).^2.*(1./r - r.^2))./2 - b.^3.*((sin(t).^2.*(1./(2.*r) - (3.*r)./2 + r.^2))./16 + (sin(t).^2.*((29.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 - 9./(128.*(g.^(1./2) + 1).^2) - 68./(512.*g.^(1./2) + 512) + 35./264))./r - r.^2.*sin(t).^2.*((11.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 + 9./(640.*(g.^(1./2) + 1).^2) - 124./(5.*(512.*g.^(1./2) + 512)) + 209./3360) - (sin(t).^4.*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./128 - 27./(256.*(g.^(1./2) + 1).^2) - 54./(512.*g.^(1./2) + 512) + 761./5632))./r + (sin(t).^2.*(1848.*r.^6.*sin(t).^2 - 2079.*r.^3.*sin(t).^2 - 7.*sin(t).^2 + 924.*r.^3 + 10))./(19712.*r.^5) - (sin(t).^2.*(5.*cos(t).^2 - 1).*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./640 - 27./(1280.*(g.^(1./2) + 1).^2) - 54./(5.*(512.*g.^(1./2) + 512)) + 829./28160))./r.^3) + b.^2.*cos(t).*sin(t).^2.*(1./r.^2 - 1);
        apply(@contour, s, 'Psi', psi, 10, +1, s.Vinf * linspace(-1, 1, 50));
        sphere();
        hold off;
        pause(1);
    end
    %}

end

function h = graph(x, y, s) 
    h = loglog(x(y>=0), y(y>=0), ['b' s], x(y<0), -y(y<0), ['r' s]);
end

function sphere()
    t = linspace(0, 2*pi, 100);
    fill(cos(t), sin(t), 'k');
end

function apply(func, sol, name, val, Rmax, dir, varargin)
    g = sol.grid.(name);
    if isempty(val)
        Z = (sol.(name));
    else
        Z = val;
    end
    I = g.x < Rmax;
    X = g.X .* cos(g.Y);
    Y = g.X .* sin(g.Y);
    X = X(I, :);
    Y = Y(I, :) * dir;
    Z = Z(I, :);
    func(X, Y, Z, varargin{:});
end


%{    
    for k = 1:numel(S.solutions)
        sol = S.solutions{k};
        s = [['\beta = ' sprintf('%.3f', sol.beta)] ...
            sprintf('  ') ...
            ['V = '  sprintf('%.3f', sol.Vinf)]];
        clf; 
        subplot(2,2,[1 3])
        colormap(jet)
        loglog(betas, V, '-b', betas, Va, ':r', ...
            sol.beta, sol.Vinf, '.b', ...
            sol.beta, sol.Vinf_approx, '.r', ...
            sol.beta, brenner(sol), '.k', ...
            'MarkerSize', 20, 'LineWidth', 2.5)
        title('Electrokinetic Solver Results', 'FontSize', 16)
        xlabel('Electric field magnitude (\beta)', 'FontSize', 12)
        ylabel('Particle velocity magnitude (V)', 'FontSize', 12)
        legend({'Numerical Solver', 'Analytical Approximation'}, ...
            'Location', 'SouthEast', 'FontSize', 8, 'Box', 'on')
        axis([0.1 10 0.001 10])
        text(0.04, 0.4, s, 'FontSize', 12);
        subplot(2,2,2);
        hold on
        show1(@contour, sol, 'Psi', 10, +1, sol.Vinf * linspace(-1, 1, 50));
        show1(@contour, sol, 'Psi', 10, -1, sol.Vinf * linspace(-1, 1, 50));
        sphere();
        hold off; 
        axis equal; axis([-4 4 -3 3]); 
        title('Fluid streamlines', 'FontSize', 14)

        subplot(2,2,4);
        hold on
        for d = [-1 1]
            show1(@surf, sol, 'C', 10, d, 'EdgeColor', 'none');
        end
        set(gca, 'CLim', [0 2])
        shading('interp')
        colorbar;
        view(2)
        sphere();
        hold off; 
        axis equal; axis([-4 4 -3 3]); 
        title('Salt Concentration', 'FontSize', 14)
        drawnow;
        print('-dpng', '-r95', sprintf('results./%03d', k-1))
    end
%}
