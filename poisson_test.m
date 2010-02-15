% function poisson_test
% Simlate and solve Poisson
Nx = 33;
Ny = 1;
sz = [Nx Ny];
u = @(x,y) exp(x);
L = @(x,y) exp(x);

x = linspace(-1, 1, Nx);
y = linspace( 0, 0, Ny);
h = x(2)-x(1);

[X, Y] = ndgrid(x, y);
U = u(X, Y);
F = L(X, Y);
F = F(2:end-1, 1:end);

% Take an initial guess.

V = U;
% V(2:end-1, 1:end) = 0;

s = [1; 0; 1];
T = 1000;

for t = 1:T
    if mod(t, 10) == 0 
        plot(x, V - u(X, Y))
        title(sprintf('%d / %d', t, T))
        pause(0.05);
    end
    V(2:end-1, 1:end) = 0.5 * ...
        (conv2(V, s, 'valid') - h^2 * F);
end
