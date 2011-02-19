function [e] = show(S)
%% Show solver's results
U0 = S.U(S.X, S.Y); % The original continuous solution
U = U0;
U(S.I) = S.Uf;
clf;

subplot 131; 
if ~isempty(S.iter_type)
    s = sprintf('(%d iterations)', S.iters);
else
    s = '(Direct)';
end
visualize(sprintf(['Solution ' s]), S.X, S.Y, U, U0);

subplot 132; 
E = U0 - U;
e = norm(E(:), inf);
visualize(sprintf('Error (L_\\infty = %.3e)', e), S.X, S.Y, E)

if isfield(S, 'residuals')
    subplot 133;
    semilogy(1:numel(S.residuals), S.residuals);
    stop_iter = find(S.residuals, 1, 'last');
    title(sprintf('Residual L_2 norm\n(after %d %s iterations)', ...
        stop_iter, S.iter_type)); 
    xlabel('Iteration #');
end

function visualize(msg, X, Y, varargin)
Z = varargin{1};
sz = size(Z).';
if sz(1) > 1 && sz(2) > 1
    surf(X, Y, Z, 'EdgeColor', 'Black');  
    xlabel('X'); 
    ylabel('Y');
else
    XY = [X(:) Y(:)]; % 1D plot hack
    t = XY(:, sz > 1);
    plot_func = @plot;
    fmt = {'b-', 'r-.'};
    hold on;
    for k = 1:numel(varargin)
        Z = varargin{k};
        plot_func(t, Z(:), fmt{k}); 
    end
    hold off;
end
title(msg);
