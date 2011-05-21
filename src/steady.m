% Steady-state solver.
% Given beta, find the velocity for zero particle force.
% Example usage:
% >> s = steady([], 1, [1 2], 3);
function [sol, V, F] = steady(sol, betas, v, iters, varargin)

    conf = [{'version', 0, 'verbose', 0}, varargin];
    function f = func(b, u) % Total force for specified (b, u)
        [sol] = force(sol, b, u, conf{:});
        f = sol.force.total;
    end

    func(betas, v(1)); % Converge near v(1) using continuation in betas.
    b = betas(end); % Last beta, for steady-state velocity.
    [~, V, F] = secant(@(u) func(b, u), v, iters, ...
        @(u, f) fprintf('B = %e, V = %e, F = %e\n', b, u, f));
end
