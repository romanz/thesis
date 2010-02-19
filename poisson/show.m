function show(X, Y, Vf, V0)
%% Show solver's results
sz = size(Vf);
subplot 121; 
if sz(1) > 1 && sz(2) > 1
    mesh(X, Y, Vf); 
    xlabel('X'); 
    ylabel('Y'); 
else
    plot([X(:) Y(:)] * (sz.' > 1), Vf(:)); % 1D plot hack
end
title('Solution'); 

subplot 122; 
E = V0 - Vf;
if sz(1) > 1 && sz(2) > 1
    mesh(X, Y, E);  xlabel('X'); ylabel('Y');
else
    plot([X(:) Y(:)] * (sz.' > 1), E(:)); % 1D plot hack
end
E = abs(E(:));
title(sprintf('Error (L_\\infty = %.2e)', max(E)))
