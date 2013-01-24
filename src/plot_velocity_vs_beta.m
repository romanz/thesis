function plot_velocity_vs_beta
    clf;
    M = {};
    hold on;
    iter(128, 10,  '-s', [1 0 0]);
    iter(128, 30,  '-o', [1 0 0]);
    [b1, v1] = iter(128, 100, '-d', [1 0 0]);
    iter(256, 10,  '-s', [0 .5 0]);
    iter(256, 30,  '-o', [0 .5 0]);
    [b2, v2] = iter(256, 100, '-d', [0 .5 0]);
    iter(512, 10,  '-s', [0 0 1]);
    iter(512, 30,  '-o', [0 0 1]);
    [b3, v3] = iter(512, 100, '-d', [0 0 1]);
    xlim([0 6])
    legend(M, 'Location', 'NorthWest')
    print -depsc2 LargeBetaV.eps
    
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
        if nargin < 3
            linespec = '-';
        end
        f = sprintf('MATs/*[%dx%d]_Rmax=%.1f*_lite.mat', N+1, N+1, Rmax);
        L = glob(f).';
        s = cellfun(@load_data, L, 'UniformOutput', false);
        res = cat(1, s{:});
        res = sortrows(res);
        beta = res(:, 1);
        V = res(:, 2);
        plot(beta, V, linespec, 'Color', c)
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