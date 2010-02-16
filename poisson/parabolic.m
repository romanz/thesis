%% Perform parabolic interpolation 
syms x0 x1 x2 A B C
y0 = A*x0^2 + B*x0 + C;
y1 = A*x1^2 + B*x1 + C;
y2 = A*x2^2 + B*x2 + C;

% a = ((y2 - y1)/(x2 - x1) - (y1 - y0)/(x1 - x0))/(x2 - x0);
a = [(x2-x1) (x0-x2) (x1-x0)]/((x2-x1)*(x1-x0)*(x2-x0)) * [y0 y1 y2].';
disp(simple(a))

p = [0 -1 1; 1 0 -1; -1 1 0] * [x0; x1; x2];
p = -p / prod(p);
p = [y0 y1 y2] * p;
disp(simple(p))
