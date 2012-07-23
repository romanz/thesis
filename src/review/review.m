clear
load clean\[400x79]_1e+07_alpha=0.00.mat
% load clean\[400x79]_1e+07_alpha=0.50.mat
% load 20120701150920.mat
r = g.Phi.r;
rD = (r(2:end)+r(1:end-1))*0.5;
phiN = regrid(sol.Phi);
% phiA_1 = -r * sol.beta;
% phiA_2 = phiA_1 + sol.beta^2 * 0.204292./r.^3;
% phiA_3 = phiA_2 + sol.beta^3 * (0.051073./r.^5-0.0276818./r.^4+2.08167*10^-17./r.^3+0.292141./r.^2);

phiA_1 = -r * sol.beta;
phiA_2 = phiA_1 + sol.beta^2 * 0.133333./r.^3;
phiA_3 = phiA_2 + sol.beta^3 * (1./(30*r.^5)+0.0171429./r.^4-0.0505556./r.^2);

semilogx(rD, diff(phiA_1)./diff(r), '-', ...
    rD, diff(phiA_2)./diff(r), '-', ...
    rD, diff(phiA_3)./diff(r), '-', ...
    rD, diff(phiN(:, 1))./diff(r), '.')

dphiN_dr = diff(phiN(:, 1))./diff(r);
dphiA_dr = diff(phiA_3)./diff(r);

legend('Linear', 'Quadratic', 'Cubic', 'Numeric')