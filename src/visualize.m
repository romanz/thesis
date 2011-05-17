%% Grab all streamline plots
load
clf;
for k = 1:numel(sols)
    sol = sols{k};
    streamlines(sol, linspace(-sol.Vinf, sol.Vinf, 100));
    axis equal;
    d = 5;
    xlim([-d d]);
    ylim([0 d])
    text(-4.5, 4.5, sprintf('Streamlines - \\beta = %f', betas(k)))
    drawnow;
    print('-dpng', sprintf('%02d', k))
end
