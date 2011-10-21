function s = solver(eqns, vars, weights, func)
    E = [];
    for w = weights
        E = [E, func(w * eqns)];
    end
    V = num2cell(vars, 1);
    s = solve(E, V{:});
end
