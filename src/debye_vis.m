clear; clc;
%%
figure(1)
r = linspace(1, 5, 201);
t = linspace(0, 2*pi, 201);
[R, T] = ndgrid(r, t);
X = R .* cos(T);
Y = R .* sin(T);
C = 1 + 0.1 * cos(T) ./ R.^2;
clf;
hold on
fill(cos(t), sin(t), 'k');
surf(X, Y, C, 'EdgeColor', 'none')
view(2)
axis equal
axis([-1 1 -1 1]*3)
m = jet;
colormap(m)
colorbar
title('Non-uniform salt concentration C = c_+ = c_-')
hold off
%%
figure(2)
Cout = 1;
gamma = exp(1); % Cp
zeta = log(Cout/gamma);
z = linspace(0, 6, 100);
phi = 4*atanh( exp(-z*sqrt(Cout))*tanh(zeta/4) );
Cp = Cout * exp(-phi);
Cm = Cout * exp(+phi);
plot(z, Cp, z, Cm, z, phi, ':');
legend({'c_+', 'c_-', '\phi'})
xlabel('\rho = (r - 1)/\delta')
text(0.2, 2.7, 'c_+ = \gamma')
text(0.2, 0.3, 'c_- = 1/\gamma')
text(5.0, 1.2, 'c_+ = c_- = C')
text(0.4, -0.8, '\zeta = \Delta \phi = log(C/\gamma)')
title('Debye-scale ionic concentrations')
