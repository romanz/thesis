% Dukhin-Derjaguin slip velocity
function V = ddslip(Phi, C, gamma, theta)
    C = col(mean(C(1:2, 2:end-1))); % on R=1
    Phi = col(mean(Phi(1:2, 2:end-1))); % on R=1
    
    h = [1;1]/2;
    xi = log(average(C(:), h)/gamma);
    lnC = log(C);
    
    D = [1;-1];
    dtheta = average(theta(:), D);
    deriv = @(f) average(f(:), D) ./ dtheta(2:end-1);
    V = xi .* deriv(Phi) + 2 * log(1 - tanh(xi/4).^2) .* deriv(lnC);
end
