function debug_solve(f, S, indices)
    syms u t real
    if ~isempty(S)
        f = subs(f, S{:});
    end
    
    f = f(:)
    y = [];
    for k = indices
        g = f;
        if k > 0
            g = diff(g, u, k);
        end
        g = subs(g, u, 0);
        y = [y; g];
    end
    y
    sol = solve(y);
    names = fieldnames(sol);
    for k = 1:numel(names)
        n = names{k};
        v = sol.(n);
        fprintf('%s = %s; ', n, char(v))
    end
    fprintf('\n')
end