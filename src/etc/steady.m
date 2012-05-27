function [sol, V, F] = steady(sol, v, iters)

    logger('steady', 'alpha = %f, beta = %e', sol.alpha, sol.beta);
    function f = solve_force(u) % Total force for specified (beta, Vinf)
        [sol] = force(sol, sol.beta, u);
        f = sol.force.total;
    end

    secant_step = secant(@(u) solve_force(u), v);
    for k = 1:iters
        [V(k), F(k)] = secant_step();
        logger('steady', '%3d) V = %e \t F = %e', k, V(k), F(k));
    end
    
end
