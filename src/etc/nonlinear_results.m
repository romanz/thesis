clear
load matlab
clf;
hold on;
bnd = @(x) mean(x(1:2, :));
S0 = solutions{1};
b = [];
v = [];
for k = 1:numel(solutions)
    S = solutions{k};
    b(k) = S.beta;
    v(k) = mean(bnd(S.Phi) / S.beta);
    plot(bnd(S.Phi) / S.beta)
    pause(0.1)
end
%% 
clf;
plot(b, v, '.-')