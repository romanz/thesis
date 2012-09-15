function plot_results_for_paper

    load data

    b = b{1};
    g = cat(2, g{:}); g = g(1, :); g = g(:);
    v = cat(2, v{:});

    I = b < 5;
    b = b(I);
    v = v(I, :);

    [~, I] = sort(abs(1-g));
    g = g(I);
    v = v(:, I);

    b_ = logspace(log10(min(b)), log10(max(b)), 1000);
    for k = 1:numel(g)
        u(:, k) = analytic(b_, g(k));
    end

    close all;
    clr = (0:3)'*0.25*[1 1 1];

    I = g <= 1;
    figure; hold on;
    set(gca, 'ColorOrder', clr)
    plot(b_, u(:, I), '-', 'LineWidth', 2);
    plot(b, v(:, I), 'o', 'MarkerSize', 8);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([min(b) max(b)])
    xlabel('Electric field magnitude (\beta)', 'FontSize', 20)
    ylabel('Steady-state velocity magnitude', 'FontSize', 20)
    L = cellfun(@(x) sprintf('\\gamma = %.4f', x), num2cell(g(I)), 'UniformOutput', false);
    legend(L, 'Location', 'SouthEast', 'FontSize', 18); legend boxoff;
    grid on;
    print('-depsc2', 'ionex_cubic1.eps')

    I = g > 1;
    figure; hold on;
    set(gca, 'ColorOrder', clr)
    plot(b_, u(:, I), '-', 'LineWidth', 2);
    plot(b, v(:, I), 'o', 'MarkerSize', 8);
    plot(b_, -u(:, I), '--', 'LineWidth', 2);
    plot(b, -v(:, I), 's', 'MarkerSize', 8);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([min(b) max(b)])
    xlabel('Electric field magnitude (\beta)', 'FontSize', 20)
    ylabel('Steady-state velocity magnitude', 'FontSize', 20)
    L = cellfun(@(x) sprintf('\\gamma = %.4f', x), num2cell(g(I)), 'UniformOutput', false);
    legend(L, 'Location', 'SouthEast', 'FontSize', 18); legend boxoff;
    grid on;
    print('-depsc2', 'ionex_cubic2.eps')
end

function v = analytic(b, g)
    B3 = (31./(320*(g.^(1./2) + 1)) - 9./(320*(g.^(1./2) + 1).^2) + 1./1680);
    B1 = 2*log(1./(2*g.^(1./2)) + 1./2);
    w3 = b.^3*(B3 - B1*11./320);
    w1 = b*B1;
    v = w1 + w3;
end