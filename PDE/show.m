function [e] = show(mat_file)
%% Show solver's results
load(mat_file)
U0 = U(X, Y); % The original continuous solution
U = U0;
U(I) = Uf;

subplot 131; 
visualize(X, Y, U, sprintf('Solution (%d iterations)', iters));

subplot 132; 
E = U0 - U;
e = norm(E(:), inf);
visualize(X, Y, E, sprintf('Error (L_\\infty = %.3e)', e))

subplot 133; 
semilogy(1:numel(residuals), residuals);
stop_iter = find(residuals, 1, 'last');
title(sprintf('Residual L_2 norm\n(after %d %s iterations)', ...
    stop_iter, iter_type)); 
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
