clc; 

for r = [10 30 100 300]
    p = sprintf('better_R_res/%d/', r);
    D = dir([p '*x*.mat']);

    if isempty(D)
        fprintf('No MAT files found.\n')
        return
    end

    S = {};
    for k = 1:numel(D)
        d = D(k);
        f = d.name;
        S{k} = load([p, f], 'betas', 'k', 'v', 'Nr');
        disp([p f])
    end

    k = numel(S{1}.betas); % number of betas
    S = reshape(S, k, numel(S)/k);
    % S = S(:, 2:end);

    betas = cellfun(@(s) s.betas(s.k), S);
    [~, I] = sort(betas(:, 1), 1);

    N = cellfun(@(s) s.Nr-1, S);
    [~, J] = sort(N(1, :), 2);

    V = cellfun(@(s) s.v(end), S);

    Vn = V(I, J);
    Vr0 = Vn(:, 1:2) * [-1/3; 4/3];
    Vr1 = Vn(:, 2:3) * [-1/3; 4/3];
    Vr2 = Vr1 + (Vr1 - Vr0)/15;

    U1 = 4.25752957407993;
    U3 = 0.966869;

    b = betas(:, 1);

    V1 = U1*b;
    V3 = U3*b.^3;

    figure(1); plot(b, [Vn(:, :)], 'o-', b, [V1+V3], 's-', b, [V1], '-')
%     legend('Numerical', 'Analytical', 'Linear', ...
%         'Location', 'SouthEast')
    grid on
    xlabel('\beta');
    ylabel('steady-state velocity: deviation from linear')
    xlim([0 6])
    ylim([0 25])
    title('Numerical results for cubic correction')
    print('-depsc2', [p 'graph.eps'])
    save([p 'data.mat'])
end

%{
figure(2); loglog(b, Vr1, 'o-', b, V1+V3, 's-', b, V1, 'd-')
legend('Numerical', 'Cubic', 'Linear', 'Location', 'SouthEast')
grid on
xlabel('\beta'); ylabel('steady-state velocity')
%}