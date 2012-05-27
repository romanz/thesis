function stationary_case

    sol = struct( ...
        'radius', logspace(0, 7, 400), ...
        'theta', linspace(0, pi, 50), ...
        'gamma', 1, ...
        'alpha', 0, ...
        'maxres', 1e-11, ...
        'iters', [2 5] ...    
    );

    function f = func(b, u)
        sol.beta = b;
        sol = force(sol, b, u);
        f = sol.force.total;
    end

    sol = force(sol); % init
    tic;

    betas = [];
    gammas = [];
    for e = 10.^(-6:0.5:-2);

        g = 1 + e;
        sol.gamma = g;

        B3 = (31/(320*(g^(1/2) + 1)) - 9/(320*(g^(1/2) + 1)^2) + 1/1680);
        B1 = 2*log(1/(2*g^(1/2)) + 1/2);
        % w = b.^3*(B3 - B1*11/320) + b*B1;
        beta_c = sqrt(-B1 / (B3 - B1*11/320));

        b = [1 2] * beta_c;
        f = [];
        step = secant(@(b) func(b, 0), b);
        for k = 1:10
            [b, f] = step()
        end
        betas = [betas; b beta_c];
        gammas = [gammas; g];
        save stat
    end

    loglog(sqrt(gammas-1), betas, '.')
end
