% [Vr, Vtheta, Psi] = stokes_solution(U, r, theta)
function [Vr, Vtheta, Psi] = stokes_solution(U, r, theta)
    Vr = U .* cos(theta).* (1 - 1.5 ./ r + 0.5 ./ (r.^3));
    Vtheta = -U .* sin(theta).* (1 - 0.75 ./ r - 0.25 ./ (r.^3));
    Psi = 0.5 * U .* r .* sin(theta) .* (1 - 1.5 ./ r + 0.5 ./ (r.^3));
end
