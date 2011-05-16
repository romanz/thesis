%% Grab all streamline plots into a movie
load
clf;
for k = 1:numel(sols)
    sol = sols{k};
    streamlines(sol, linspace(-sol.Vinf, sol.Vinf, 100));
    axis equal;
    d = 5;
    xlim([-d d]);
    ylim([0 d])
    title(sprintf('Streamlines - \\beta = %f', betas(k)))
    drawnow;
    F(k) = getframe;    
end
%%
movie(F)
