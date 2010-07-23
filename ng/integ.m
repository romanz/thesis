function S = integ(f, a, b, n)
h = (b - a) / (2*n); % interval length
x = linspace(a, b, 2*n + 1); % Newton-Cotes nodes
F = f(x); % Function values
I = logical(mod(1:2*n-1, 2)); % [1 0 1 0 ... 0 1] pattern
I = [1, (mod(1:2*n-1, 2) + 1) * 2, 1]; % [1 4 2 4 2 ... 2 4 1] pattern
w = I * h / 6; % Simpson's weights
S = w * F(:);
