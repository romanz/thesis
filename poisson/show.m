function [e] = show(mat_file)
%% Show solver's results
load(mat_file)
V = V0;
V(I) = Vf;

subplot 131; 
visualize(X, Y, V, sprintf('Solution (%d iterations)', iters));

subplot 132; 
E = V0 - V;
e = norm(E(:), inf);
visualize(X, Y, E, sprintf('Error (L_\\infty = %.3e)', e))

subplot 133; 
semilogy(1:iters, residuals);
stop_iter = find(residuals, 1, 'last');
title(sprintf('Residual L_2 norm\n(after %d iterations)', stop_iter)); 
xlabel('Iteration #');

function visualize(X, Y, Z, t)
sz = size(Z).';
if sz(1) > 1 && sz(2) > 1
    surf(X, Y, Z, 'EdgeColor', 'Black');  
    xlabel('X'); 
    ylabel('Y');
else
    XY = [X(:) Y(:)]; % 1D plot hack
    plot(XY(:, sz > 1), Z(:), '-'); 
end
title(t);
