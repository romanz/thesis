function plot_velocity_vs_beta
    clf;
    M = {};
    hold on;
    iter(128, 10,  '-s', [1 0 0])
    iter(128, 30,  '-o', [1 0 0])
    iter(128, 100, '-d', [1 0 0])
    iter(256, 10,  '-s', [0 .5 0])
    iter(256, 30,  '-o', [0 .5 0])
    iter(256, 100, '-d', [0 .5 0])
    iter(512, 10,  '-s', [0 0 1])
    iter(512, 30,  '-o', [0 0 1])
    iter(512, 100, '-d', [0 0 1])
    xlim([0 6])
    legend(M, 'Location', 'NorthWest')
    print -depsc2 LargeBetaV.eps

    function iter(N, Rmax, linespec, c)
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