function vorticity_test
force_solver('results', 'w', 0, 2.2, 1e-2, 1e3, [100 30], 30, [false 100 false]);
load results
[W, gW] = vorticity(solVx, solVy, gridVx, gridVy);
Wt = 1.5 * Vinf * sin(gW.Y) ./ gW.X.^2;
% clf;
% contour(gW.y, gW.x, W, 20)
% ylim([1 10])
subplot 211; mesh(W); zlim([0 1.5]*Vinf); title('Numeric results')
subplot 212; mesh(Wt); zlim([0 1.5]*Vinf); title('Theory')
E = W - Wt; 
fprintf('Vorticity error: %.2f%%\n', 100 * norm(E(:), inf) / norm(W(:), inf))