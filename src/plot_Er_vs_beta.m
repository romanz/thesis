function plot_Er_vs_beta
    clf;
    M = {};
    hold on;
    for b = 0.1:0.1:10
        iter(b, 128,  10, '-d', [0 0 1]);
        iter(b, 128,  30, '-d', [0 0.5 0]);
        iter(b, 128, 100, '-d', [1 0 0]);
        iter(b, 256,  10, '-s', [0 0 1]);
        iter(b, 256,  30, '-s', [0 0.5 0]);
        iter(b, 256, 100, '-s', [1 0 0]);
        iter(b, 512,  10, '-o', [0 0 1]);
        iter(b, 512,  30, '-o', [0 0.5 0]);
        iter(b, 512, 100, '-o', [1 0 0]);
        drawnow;
    end
%     legend(M, 'Location', 'NorthWest')
    print -depsc2 LargeBetaV.eps

    function iter(b, N, Rmax, linespec, c)
        0;
        if nargin < 4
            linespec = '-';
        end
        f = sprintf('MATs/*beta=%.3e_[%dx%d]_Rmax=%.1f*0.mat', b, N+1, N+1, Rmax);
        L = glob(f).';
        if isempty(L)
            1;
            return;
        end
        f = L{1};
        s = load(f);
        Phi = s.sol.Phi;
        r = Phi.grid.r;
        Phi = regrid(Phi);
        Er = -diff(Phi(:, 1)) ./ diff(r);
        r = convn(r, [1;1]/2, 'valid');
        plot(b, Er(1) / b, linespec, 'Color', c)
        M{numel(M)+1} = sprintf('[%dx%d] R_{max}=%.0f', N, N, Rmax);
%         xlim([1 2])
    end

end

function res = load_data(f)
    s = load(f);
    s = s.sol;
    res = [s.beta, s.Vinf];
end