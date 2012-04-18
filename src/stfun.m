function r = stfun(s, names, f)
    r = struct();
    for k = 1:numel(names)
        n = names{k};
        v = s.(n);
        r.(n) = f(v);
    end
end
