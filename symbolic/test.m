clear;
load assert
syms r t b A B pi real
[num, den] = numden(simplify(e / b^3));
s = expand(num);
pretty(s)

for k = 1:numel(s)
    
end

% [C, P] = coeffs(s, r);
% 
dot = @(f, g) int(f*g, t, 0, pi) / pi;
S = {dot(s(1), cos(t)); dot(s(1), cos(t)^3);
     dot(s(2), cos(t)); dot(s(2), cos(t)^3)};
sol = solve(S{:}, 'A1', 'A2', 'A3', 'A4')

% power = @(f, r) subs(diff(f, r) / f, r, 1);
% 
% sol = 0;
% for k = 1:numel(C)
%     c = C(k);
%     p = P(k);
%     u = (A*cos(t) + B*cos(t)^3) * p * r^2 / den;
%     Lu = scalar_laplacian(u);
%     
%     f = c * p / den;
%     res = simple((f - Lu) * den / p);
%     eq(1) = dot(res, cos(t));
%     eq(2) = dot(res, cos(t)^3);
%     
%     res = subs(eq(:), {A,B}, {0,0});
%     H = [ diff(eq(:),A) diff(eq(:),B) ];
%     
%     I = [ diff(Lu, A); diff(Lu,B)] ~= 0;
%     R = eye(2); 
%     R = R(I, :);
%     
%     s = R'*(R*H*R' \ R*res);
%     sol = sol + subs(u, [A,B], s);
% end
% sol = simple(sol);
% pretty(sol)