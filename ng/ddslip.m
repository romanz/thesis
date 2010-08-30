% Dukhin-Derjaguin slip velocity
function V = ddslip(solPhi, solC, gamma, theta)
    C = col(mean(solC(1:2, 2:end-1))); % on R=1
    Phi = col(mean(solPhi(1:2, 2:end-1))); % on R=1
    
    xi = log(average(C(:), [1;1]/2)) - log(gamma);
    lnC = log(C);
    
    dtheta = diff(theta(:));
    deriv = @(f) diff(f(:)) ./ dtheta;
    V = xi .* deriv(Phi) + 2 * log(1 - tanh(xi/4).^2) .* deriv(lnC);
end
