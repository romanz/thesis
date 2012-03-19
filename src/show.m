function show(name, first)
    if nargin < 2
        first = 1;
    end
    load(name);
    figure(1); clf;
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'Position',[100 100 600 800]);
    try
        v = attr(solutions, 'Vinf');
        betas = attr(solutions, 'beta');
    catch err
        disp(err)
        return
    end
    disp(name)

    g = sol.gamma;
    B3 = (31./(320*(g.^(1./2) + 1)) - 9./(320*(g.^(1./2) + 1).^2) + 1./1680);
    B1 = 2*log(1./(2*g.^(1./2)) + 1./2);
    b = logspace(log10(min(betas)), log10(max(betas)), 1000);
    a = sol.alpha;
    w = (248./(2560.*g.^(1./2) + 2560) - 9./(320.*(g.^(1./2) + 1).^2) - (a.*log(1./(2.*g.^(1./2)) + 1./2))./8 - (11.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./160 - (61.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./80 - (83.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./160 + (176.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(1280.*g.^(1./2) + 1280) + 1./1680).*b.^3 + 2.*log(1./(2.*g.^(1./2)) + 1./2).*b;
    w3 = b.^3*(B3 - B1*11./320);
    w1 = b*B1;
    clf; 
    hold on; 
    graph(betas, v, 'o', 's'); 
    graph(b, w1+w3, '-', '-.'); 
%     graph(b, w1, '-', '-'); 
%     graph(b, w3, '-', '-'); 
    filename = sprintf('grid = [%d x %d], Rmax = %.1e, gamma = %.4f, alpha = %.4f', ...
        numel(sol.grid.radius), numel(sol.grid.theta), max(sol.grid.radius), ...
        sol.gamma, sol.alpha);
    filename = sprintf('g=%f', g);
    title(sprintf('\\gamma = %f', g), 'FontSize', 16, 'FontWeight', 'bold')
    disp(filename)
    xlabel('\beta (electric field)', 'FontSize', 16)
    ylabel('U (steady-state velocity)', 'FontSize', 16)
    hold off;
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'FontWeight', 'bold')
    ylim([1e-7 1e2])
    xlim([1e-2 5])
    grid on; grid minor;
    s = [filename '.eps'];
    print('-deps2', s); disp(s)
    if first, return; end
    betas = betas(betas <= 0.1);
    for k = 1:numel(betas)
        s = solutions{k};
        b = s.beta;
        r = s.grid.Psi.X;
        t = s.grid.Psi.Y;
        g = s.gamma;
        
        figure(2); clf; hold on        
        psi = b.^3.*(r.^2.*sin(t).^2.*((11.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 + 9./(640.*(g.^(1./2) + 1).^2) - 124./(5.*(512.*g.^(1./2) + 512)) + (61.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./160 + (83.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./320 - (88.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(5.*(256.*g.^(1./2) + 256)) + 209./3360) - (sin(t).^2.*(1848.*r.^6.*sin(t).^2 - 2079.*r.^3.*sin(t).^2 - 7.*sin(t).^2 + 924.*r.^3 + 10))./(19712.*r.^5) - (sin(t).^2.*((29.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 - 9./(128.*(g.^(1./2) + 1).^2) - 68./(512.*g.^(1./2) + 512) + (59.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./80 + (11.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./20 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./40 - (104.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(256.*g.^(1./2) + 256) + 35./264))./r + (sin(t).^4.*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./128 - 27./(256.*(g.^(1./2) + 1).^2) - 54./(512.*g.^(1./2) + 512) + (57.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./128 + (93.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./256 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./32 - (108.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(256.*g.^(1./2) + 256) + 761./5632))./r + (sin(t).^2.*((a.*log(1./(2.*g.^(1./2)) + 1./2))./8 - 1./8).*(1./(2.*r) - (3.*r)./2 + r.^2))./2 + (sin(t).^2.*(5.*cos(t).^2 - 1).*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./640 - 27./(1280.*(g.^(1./2) + 1).^2) - 54./(5.*(512.*g.^(1./2) + 512)) + (57.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./640 + (93.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./1280 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./160 - (108.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(5.*(256.*g.^(1./2) + 256)) + 829./28160))./r.^3) + b.*sin(t).^2.*log(1./(2.*g.^(1./2)) + 1./2).*(1./r - r.^2) - b.^2.*cos(t).*sin(t).^2.*((3.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 9./(16.*g.^(1./2) + 16)).*(1./r.^2 - 1);
        I = true(s.grid.Psi.sz); I(:, [1 end]) = false;
        psi(I) = psi(I) ./ (r(I) .* sin(t(I)));
        apply(@contour, s, 'Psi', psi, 15, s.Vinf * linspace(-1, 1, 35));
        sphere();
        axis equal;
        axis([-1 1 -.5 .5]*10)        
        hold off;
        txt = sprintf('\\beta = %.2f', betas(k));
        title(txt, 'FontSize', 16, 'FontWeight', 'bold');
        set(gca, 'XTick', [], 'YTick', [])
        n = sprintf('Stream/A%.2f.eps', b);
        print('-deps2', n)

        fprintf('\\includegraphics[width=0.45\\textwidth]{figs/%s}\n', n)

        figure(3); clf; hold on        
        psi = [];
        apply(@contour, s, 'Psi', psi, 15, s.Vinf * linspace(-1, 1, 35));
        sphere();
        axis equal;
        axis([-1 1 -.5 .5]*10)        
        hold off;
        txt = sprintf('\\beta = %.2f', betas(k));
        title(txt, 'FontSize', 16, 'FontWeight', 'bold');
        set(gca, 'XTick', [], 'YTick', [])
        n = sprintf('Stream/N%.2f.eps', b);
        print('-deps2', n)
        
        fprintf('\\includegraphics[width=0.45\\textwidth]{figs/%s}\n', n)
    end
    %}

end

function graph(x, y, pos, neg) 
    plot(x(y>=0), y(y>=0), pos, 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 2)
    plot(x(y<0), -y(y<0), neg, 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 2);
end

function sphere()
    t = linspace(0, 2*pi, 100);
    fill(cos(t), sin(t), 'k');
end

function apply(func, sol, name, val, Rmax, varargin)
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
    Y = Y(I, :);
    Z = Z(I, :);
    func([X, X], [Y, -Y], [Z, -Z], varargin{:});
%     func([X, X], [Y, -Y], [Z, -Z], varargin{:});
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
