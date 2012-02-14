clc
clear;
load assert
syms r t b pi real
e = e(1:2);
[num, den] = numden(e);
num = expand(num);
C = [];
for k = 1:numel(num)
    c = num(k);
    [c, p] = coeffs(c, r);
    c = c.*p / den(k);
    C = [C, (c).'];
end
sol = 0;
for k = 1:size(C ,1)
    f = C(k, :).';
    p = (r * diff(f(1), r) ./ f); % powers of r
    p = r ^ p(1);

    L = @(psi) simple(vector_laplacian(curl(psi)));
    G = @(p) -simple(grad(p));

    ops = {L, L, G, G};
    u = [sin(t)^2*[1, sin(t)^2]*r^4; cos(t)*r*[1, sin(t)^2]].'*p;

    dot = @(f, g) int(f .* g, t, 0, pi) / pi;

    v = [];
    for j = 1:numel(u)
        op = ops{j};
        v = [v, op(u(j))];
    end
    
    f = expand(b^-3 * (f/p));
    v = (v/p);
    
    w = [cos(t) cos(3*t); sin(t) sin(3*t)];
    M = [];
    
    for i = 1:size(w, 1)
        for j = 1:size(w, 2)
            M = [M; dot([v(i, :) f(i)], w(i, j))];
        end
    end
    R = rref(M)
    I = (R(:, end) ~= 0);
    [~, J] = find(R(I, 1:end-1) ~= 0);
    
    c = sym(zeros(numel(u), 1));
    c(J) = R(I, end);
    sol = sol + sum(u .* reshape(c, size(u)));
end

sol = simple(sol.')