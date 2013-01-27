function plot_Er_layer
clear
load results_Er(beta)_512x512_Rmax=10.mat
b = betas(1:end-1)
F = (E - repmat(E(end, :), size(E, 1), 1));
F = F * diag(1./max(F));
for j = 1:size(F, 2)
    f = F(:, j);
    q(j) = interp1(f, r, 0.5);
end
clf;
subplot(1,2,1)
hold on;
plot(r, F, '-');
plot(q, 0.5, '.');
xlim([1 1.2])

subplot(1,2,2)
hold on
y = log10(q-1);
I = 10:60;
r = linreg(b(I), y(I));
plot(b, y, '.', b, r.a*b + r.b)
end