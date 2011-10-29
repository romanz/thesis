function f = series(f, x, x0, n)
    f = simple(taylor(f, n+1, x, x0));
end
