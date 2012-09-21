function f = series(f, x, x0, n)
    f = simple(f);
    if any(f ~= 0)  
        t = taylor(f, x, 'ExpansionPoint', x0, 'Order', n+1);
        f = simple(t);
    end
end
