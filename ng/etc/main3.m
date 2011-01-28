figure(1);
clf;
subplot 121
hold on
for k = [-6:-1]
    f = sprintf('gamma%d', k);
    s = load(f);
    loglog(log10(s.betas), log10(s.Vinf), log10(s.betas), log10(s.U), '.')
end
subplot 122
hold on;
for k = [1:6]
    f = sprintf('gamma%d', k);
    s = load(f);
    plot(log10(s.betas), log10(-s.Vinf), log10(s.betas), log10(-s.U), '.')
end
hold off;