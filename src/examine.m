function examine
clf;
show('20120715154948', '--', 30); % 150 o
show('20120715153336', ':', 20); % 100 :
show('20120715153830', '-', 10); % 50 -
plot([0 30], [0 0])

function show(fname, marker, k_index)
load(fname)
hold on

sol.beta * (sol.Du*log(16) + sol.zeta)/(1 + 2*sol.Du) - sol.Vinf

r = g.Phi.r;
rD = (r(2:end)+r(1:end-1))*0.5;
phiN = regrid(sol.Phi);

mu = cos(g.Phi.t(k_index));
phiA_1 = -r * mu .* sol.beta;
phiA_2 = phiA_1 + sol.beta^2 * ((0. +(0.0666667*(-1+3*mu.^2))./r.^3+(0.166667 + (1/12)*(1-3*mu.^2))./r));
phiA_3 = phiA_2 + sol.beta^3 * ((0. +(0.00857143*(-3*mu+5*mu^3))./r.^4+(mu/20+(1/40)*(3*mu-5*mu^3))./r.^3+(-0.0605556*mu+(1/200)*(-3*mu+5*mu^3))./r.^2+(mu/75+(1/100)*(-3*mu+5*mu^3))./r.^5));

% semilogx(rD, diff(phiA_1)./diff(r), '-', ...
%     rD, diff(phiA_2)./diff(r), '-', ...
%     rD, diff(phiA_3)./diff(r), '-', ...
%     rD, diff(phiN(:, k))./diff(r), '.')

semilogx(rD, diff(phiA_3)./diff(r) - diff(phiN(:, k_index))./diff(r), marker)
% 
% dphiN_dr = diff(phiN(:, k))./diff(r);
% dphiA_dr = diff(phiA_3)./diff(r);
% 
% legend('Linear', 'Quadratic', 'Cubic', 'Numeric')
