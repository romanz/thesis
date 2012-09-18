function plot_stream_for_paper(name)
    S = load(name);
    fprintf('%s\n', name);
    solutions = S.solutions;

    for k = 1:numel(solutions)
        s = solutions{k};
        b = s.beta;
        r = s.grid.Psi.X;
        t = s.grid.Psi.Y;
        g = s.gamma;
        a = s.alpha;
        
        figure(1); clf; hold on        
        psi = b.^3.*(r.^2.*sin(t).^2.*((11.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 + 9./(640.*(g.^(1./2) + 1).^2) - 124./(5.*(512.*g.^(1./2) + 512)) + (61.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./160 + (83.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./320 - (88.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(5.*(256.*g.^(1./2) + 256)) + 209./3360) - (sin(t).^2.*(1848.*r.^6.*sin(t).^2 - 2079.*r.^3.*sin(t).^2 - 7.*sin(t).^2 + 924.*r.^3 + 10))./(19712.*r.^5) - (sin(t).^2.*((29.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./320 - 9./(128.*(g.^(1./2) + 1).^2) - 68./(512.*g.^(1./2) + 512) + (59.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./80 + (11.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./20 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./40 - (104.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(256.*g.^(1./2) + 256) + 35./264))./r + (sin(t).^4.*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./128 - 27./(256.*(g.^(1./2) + 1).^2) - 54./(512.*g.^(1./2) + 512) + (57.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./128 + (93.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./256 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./32 - (108.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(256.*g.^(1./2) + 256) + 761./5632))./r + (sin(t).^2.*((a.*log(1./(2.*g.^(1./2)) + 1./2))./8 - 1./8).*(1./(2.*r) - (3.*r)./2 + r.^2))./2 + (sin(t).^2.*(5.*cos(t).^2 - 1).*((9.*log((g.^(1./2) + 1)./(2.*g.^(1./2))))./640 - 27./(1280.*(g.^(1./2) + 1).^2) - 54./(5.*(512.*g.^(1./2) + 512)) + (57.*a.^2.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2).^2)./640 + (93.*a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*log(1./(2.*g.^(1./2)) + 1./2))./1280 + (a.*log((g.^(1./2) + 1)./(2.*g.^(1./2))).*((27.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 81./(16.*g.^(1./2) + 16)))./160 - (108.*a.*log(1./(2.*g.^(1./2)) + 1./2))./(5.*(256.*g.^(1./2) + 256)) + 829./28160))./r.^3) + b.*sin(t).^2.*log(1./(2.*g.^(1./2)) + 1./2).*(1./r - r.^2) - b.^2.*cos(t).*sin(t).^2.*((3.*log(1./(2.*g.^(1./2)) + 1./2).*(2.*a.*log(1./(2.*g.^(1./2)) + 1./2) + 1))./8 - 9./(16.*g.^(1./2) + 16)).*(1./r.^2 - 1);
        I = true(s.grid.Psi.sz); I(:, [1 end]) = false;
        psi(I) = psi(I) ./ (r(I) .* sin(t(I)));
        apply(@contour, s, 'Psi', psi, 15, s.Vinf * linspace(-1, 1, 35));
        sphere();
        axis equal;
        axis([-1 1 -.5 .5]*10)        
        hold off;
        txt = sprintf('Analytical: \\alpha = %.1f, \\beta = %.2f, \\gamma = %.4f', a, b, g);
        title(txt, 'FontSize', 16, 'FontWeight', 'bold');
        set(gca, 'XTick', [], 'YTick', [])
        n = sprintf('Stream/A%_4f_%.4f.eps', k, g, b);
        print('-deps2', n)

        %fprintf('\\includegraphics[width=0.45\\textwidth]{figs/%s}\n', n)

        figure(2); clf; hold on        
        psi = [];
        apply(@contour, s, 'Psi', psi, 15, s.Vinf * linspace(-1, 1, 35));
        sphere();
        axis equal;
        axis([-1 1 -.5 .5]*10)        
        hold off;
        txt = sprintf('Numerical: \\alpha = %.1f, \\beta = %.2f, \\gamma = %.4f', a, b, g);
        title(txt, 'FontSize', 16, 'FontWeight', 'bold');
        set(gca, 'XTick', [], 'YTick', [])
        n = sprintf('Stream/N%_4f_%.4f.eps', k, g, b);
        print('-deps2', n)
        
        %fprintf('\\includegraphics[width=0.45\\textwidth]{figs/%s}\n', n)
    end

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
end
