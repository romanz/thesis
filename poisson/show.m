function [e] = show(mat_file)
%% Show solver's results
load(mat_file)

subplot 121; 
visualize(X, Y, Vf, sprintf('Solution (%d iterations)', iters));

subplot 122; 
E = V0 - Vf;
e = norm(E(:), inf);
visualize(X, Y, E, sprintf('Error (L_\\infty = %.3e)', e))

function visualize(X, Y, Z, t)
sz = size(Z).';
if sz(1) > 1 && sz(2) > 1
    surf(X, Y, Z, 'EdgeColor', 'Black');  
    xlabel('X'); 
    ylabel('Y');
else
    plot([X(:) Y(:)] * (sz > 1), Z(:), '-'); % 1D plot hack
end
title(t);
