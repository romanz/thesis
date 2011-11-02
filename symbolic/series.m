function f = series(f, x, x0, n)
    f = simple(f);
    if any(f ~= 0)       
        f = simple(taylor(f, n+1, x, x0));
    end
end
