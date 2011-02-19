
radius = [1 2 4 8]';
theta = [0:3]';
[center, interior, xstag, ystag] = ...
    grids(radius, theta);
gridPhi = center;
gridVx = xstag;
gridVy = ystag;
gridP = interior;
gridC = center;

clf;

xrect = [1 8 8 1 1];
yrect = [3 3 0 0 3];

plot(gridC.X(:), gridC.Y(:), 'r.', ...
gridPhi.X(:), gridPhi.Y(:), 'd', ...
gridP.X(:), gridP.Y(:), 's', ...
gridVx.X(:), gridVx.Y(:), '>k', ...
gridVy.X(:), gridVy.Y(:), '^k', ...
xrect, yrect, '--g', 'LineWidth', 1.5)
set(gca, 'XTick', radius, 'YTick', theta)
set(gca, 'XGrid', 'on', 'YGrid', 'on')
axis([0 12 -1 4])
set(gca, 'FontSize', 14)
legend('C', '\Phi', 'P', 'V_R', 'V_\theta', 'Boundary')
xlabel('R', 'FontSize', 14); ylabel('\theta', 'FontSize', 14)

print -depsc2 StaggeredGrid.eps


% F = @(x) x.^3;
% 
% r0 = 1;
% dr = 0.1;
% r = r0 + [-dr*0.9, 0, dr];
% f = F(r);
% M = [r.^0; r.^1; r.^2];
% a = f / M
% g = @(x) a(1) + a(2)*x + a(3)*x.^2;
% 
% y = 2*r0*(a(2) + 3*a(3)*r0)
% z = (0.25*(r(2)+r(1))^2*((f(2)-f(1)) / (r(2)-r(1))) - ...
%      0.25*(r(3)+r(2))^2*((f(3)-f(2)) / (r(3)-r(2)))) / ((r(3)-r(1))/2)
% 
% x = linspace(min(r), max(r), 100);
% %plot(r, f, '.', x, F(x), '-', x, g(x), ':')
% plot(r, 0*f, '.', x, F(x) - g(x), '-')
% 
% 
% % dr = 10.^[-3:0.01:-1]
% % x = 0*dr;
% % for i = 1:numel(dr)
% %     x(i) = func(dr(i));
% % end
% % loglog(dr, abs(x / x(end)), '.-')
% % 
% % function x = func(dr)
% % r = 1;
% % f = @(x) 1./x;
% % 
% % dr1 = dr*0.6; f1 = ((f(r+dr1)-f(r))/dr1);
% % dr0 = dr; f0 = ((f(r)-f(r-dr0))/dr0);
% % x = ((r+dr1/2)^2 * f1 - (r-dr0/2)^2 * f0) / ((dr1+dr0)/2);
% % % fprintf('%e\n', x)