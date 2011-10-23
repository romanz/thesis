clear;
load assert
syms r t real
v = [];
A = [];
for k = [1 3]
    for n = [0 2 3 5 6]
        if n == 2 && k == 1, continue; end
        s = sprintf('A%d_%d', k, n);
        a = sym(s, 'real');
        v = [v; a * r^-n * cos(t)^k];
        A = [A; a];
    end
end
for i = 1:numel(v)
    u(i) = scalar_laplacian( v(i) );
end
s = e - sum(u);
s = numden(s);
C = coeffs(s, r);
E = [];
for c = C
    E = [E coeffs(c, cos(t))];
end

A = num2cell(A);
sol = solve(E, A{:});
res = simple(subs(s, A, struct2cell(sol)))