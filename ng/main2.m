function main(x)

betas = linspace(1e-4, 1e-3, 10);
gamma = exp(x);

U = zeros(size(betas));
for k = 1:numel(betas)
    Vinf = betas(k) * 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));
    Vinf = Vinf * linspace(0.9, 1.1, 3);

    N = 2e3;
    F = zeros(size(Vinf));
    for i = 1:numel(Vinf)
        F(i) = main1('results', 1, N, betas(k), gamma, Vinf(i));
    end

    [a0, b0, e, R] = linreg(F(:), Vinf(:));
    b0, R
    U(k) = b0;
    main1('results', 1, N, betas(k), gamma, U(k));
end
Vinf = betas * 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25))

f = sprintf('gamma[%d]', x);
save(f);