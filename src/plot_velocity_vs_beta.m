function plot_velocity_vs_beta
    Du = 1;
    zeta = 10;
    beta = linspace(0, 6);
    U = beta * (Du*log(16) + zeta)/(1 + 2*Du);
    
    figure(1)
    clf; hold on; M = {};
    iter(128, 10,  '-o', [1 0 0]);
    iter(256, 10,  '-o', [0 .5 0]);
    [b, v] = iter(512, 10,  '-o', [0 0 1]);
    plot(beta, U, '--k');
    M{end+1} = ['Linear Solution (small \beta): ' sprintf('%.3f \\beta', (Du*log(16) + zeta)/(1 + 2*Du))];
    r = linreg(b(10:40), v(10:40))
    plot(beta, r.a * beta + r.b, '-.k')
    xlim([0 6])
    ylim([0 30])
    M{end+1} = ['Linear Fit (moderate \beta): ' sprintf('%.3f \\beta %.3f', r.a, r.b)];
    legend(M, 'Location', 'NorthWest')
    xlabel('\beta')
    ylabel('Steady-state velocity')
    print -depsc2 LargeBetaV_[128_256_512]_Rmax=10.eps

    figure(2)
    clf; hold on; M = {};
    iter(512, 100, '-o', [1 0 0]);
    iter(512, 30,  '-o', [0 .5 0]);
    [b, v] = iter(512, 10,  '-o', [0 0 1]);
    plot(beta, U, '--k');
    M{end+1} = ['Linear Solution (small \beta): ' sprintf('%.3f \\beta', (Du*log(16) + zeta)/(1 + 2*Du))];
    r = linreg(b(10:40), v(10:40))
    plot(beta, r.a * beta + r.b, '-.k')
    xlim([0 6])
    ylim([0 30])
    M{end+1} = ['Linear Fit (moderate \beta): ' sprintf('%.3f \\beta %.3f', r.a, r.b)];
    legend(M, 'Location', 'NorthWest')
    xlabel('\beta')
    ylabel('Steady-state velocity')
    print -depsc2 LargeBetaV_[512x512]_Rmax=10_30_100.eps
    
    figure(3)
    clf; hold on; M = {};
    [b1, v1] = iter(128, 100, '-d', [1 0 0]);
    [b2, v2] = iter(256, 100, '-d', [0 .5 0]);
    [b3, v3] = iter(512, 100, '-d', [0 0 1]);
    clf;
    x_max = 5;
    [b1, v1] = select(b1, v1, @(x) x < x_max);
    [b2, v2] = select(b2, v2, @(x) x < x_max);
    [b3, v3] = select(b3, v3, @(x) x < x_max);
    plot(b1, [(v2 - v1) ./ (v3 - v2)], '.-');
    xlabel('\beta')
    ylabel('(V_{128x128} - V_{256x256}) / (V_{256x256} - V_{512x512})')
    print -depsc2 LargeBetaRatio.eps

    function [beta, V] = iter(N, Rmax, linespec, c)
        f = sprintf('MATs/*[%dx%d]_Rmax=%.1f*_lite.mat', N+1, N+1, Rmax);
        L = glob(f).';
        s = cellfun(@load_data, L, 'UniformOutput', false);
        res = cat(1, s{:});
        res = sortrows(res);
        beta = res(:, 1);
        V = res(:, 2);
        if ~isempty(c) && ~isempty(linespec)
            plot(beta, V, linespec, 'Color', c)
        end
        M{numel(M)+1} = sprintf('[%dx%d] R_{max}=%.0f', N, N, Rmax);
    end

end

function res = load_data(f)
    s = load(f);
    s = s.sol;
    res = [s.beta, s.Vinf];
end

function [x, y] = select(x, y, pred)
    I = pred(x);
    x = x(I);
    y = y(I);
end