function plot_Psi_symmetry
clf;
L = glob('MATs/Psi/*.mat');
res = cellfun(@show, L, 'UniformOutput', false);
res = cat(1, res{:});
res = sortrows(res);
plot(res(:, 1), res(:, 2), '.-')
xlabel('\beta')
ylabel('Relative streamlines asymmetry')
xlim([0 5])
ylim([0 0.01])
print -depsc2 streamlines_asymmetry

function res = show(f)
s = load(f);
err = s.z - fliplr(s.z);
err = norm(err(:)) / norm(s.z);
% clf; hold on;
% contour(s.x, s.y, s.z, s.c)
% contour(s.x, -s.y, s.z, s.c)
% t = linspace(0, 2*pi, 100);
% fill(cos(t), sin(t), [1 1 1]*0.5)
% axis([-1 1 -1 1]*3)
% axis equal
% pause(.1)
res = [s.beta, err];