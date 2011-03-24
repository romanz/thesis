function A = interpolator(Z, Zi)
    dir = size(Z) - size(Zi);
    assert(all(sort(dir) == [0 1]))

    K = true(size(Z));
    P1 = select(shift(K, -dir));
    P2 = select(shift(K, +dir));
    d = P2*Z(:) - P1*Z(:); assert(all(d(:) > 0));
    
    a1 = (P2 * Z(:) - Zi(:)) ./ d(:);
    a2 = (Zi(:) - P1 * Z(:)) ./ d(:);
    A = spdiag(a1) * P1 + spdiag(a2) * P2;
    assert(all(A(:) <= 1));
    assert(all(A(:) >= 0));
end
