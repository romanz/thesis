function show(mat_file)
%% Show solver's results
load(mat_file)

subplot 121; 
mesh_plot(X, Y, Vf, sprintf('Solution (%d iterations)', T));

subplot 122; 
E = V0 - Vf;
e = norm(E(:), inf);
mesh_plot(X, Y, E, sprintf('Error (L_\\infty = %.2e)', e))

function mesh_plot(X, Y, Z, t)
sz = size(Z).';
if sz(1) > 1 && sz(2) > 1
    mesh(X, Y, Z);  
    xlabel('X'); 
    ylabel('Y');
else
    plot([X(:) Y(:)] * (sz > 1), Z(:)); % 1D plot hack
end
title(t);
