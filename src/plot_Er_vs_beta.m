function plot_Er_vs_beta
    clf;
    M = {};
    hold on;
    index = 0;
    E = [];
    C = [];
    [~, ~, r] = iter(0.1, 512,  10, '-o', [0 0 1]);
    betas = 0.1:0.1:10;
    for b = betas
        index = index + 1;
        [e, c] = iter(b, 512,  10, '-o', [0 0 1]);
        if isempty(e) break; end
        E(:, index) = e;
        C(:, index) = c;
    end
%     legend(M, 'Location', 'NorthWest')
    betas = betas(1:index);
    save results_Er(beta)_512x512_Rmax=10.mat E r betas
    save results_C(beta)_512x512_Rmax=10.mat C r betas

    function [Er, C, r] = iter(b, N, Rmax, linespec, c)
        Er = []; r = []; C = [];
        if nargin < 4
            linespec = '-';
        end
        f = sprintf('MATs/*beta=%.3e_[%dx%d]_Rmax=%.1f*0.mat', b, N+1, N+1, Rmax);
        L = glob(f).';
        if isempty(L)
            return;
        end
        f = L{1};
        disp(f)
        s = load(f);
        Phi = s.sol.Phi;
        r = Phi.grid.r;
        Phi = regrid(Phi);
        Er = -diff(Phi(:, 1)) ./ diff(r);
        
        C = s.sol.C;
        C = regrid(C);
        C = C(:, 1);
        r = convn(r, [1;1]/2, 'valid');
        C = convn(C, [1;1]/2, 'valid');
%         plot(b, Er(1) / b, linespec, 'Color', c)
%         M{numel(M)+1} = sprintf('[%dx%d] R_{max}=%.0f', N, N, Rmax);
%         xlim([1 2])
    end

end

function res = load_data(f)
    s = load(f);
    s = s.sol;
    res = [s.beta, s.Vinf];
end