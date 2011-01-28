function test_1d
    clc;
    R = 10;
    N = 100;
    r = logspace(log10(1), log10(1+R), N);
    % r = linspace((1), (1+L), N);
    r = r(:);
    f = 1 ./ [1, R+1];
    L = lapl(r);
    A = L(:, 2:end-1);
    b = -f(1)*L(:, 1)-f(2)*L(:, end);
    x = A \ b;
    x = full([f(1); x; f(2)]);
    y = 1./r;
    plot(r, x - y, '.')
    norm(L * y, inf)
end

function L = lapl(x)
    I = interior(size(x));
    n = nnz(I);
    m = numel(I);
    
    K = find(I);
    Kl = find(shift(I, [-1 0]));
    Kr = find(shift(I, [+1 0]));
    
    wl = ((x(K) + x(Kl))/2).^2 ./ (x(K) - x(Kl));
    wr = ((x(Kr) + x(K))/2).^2 ./ (x(Kr) - x(K));
    w = 2./(x(Kr) - x(Kl));
    Dl = sparse([1:n, 1:n], [Kl K], wl*[-1 1], n, m);
    Dr = sparse([1:n, 1:n], [K Kr], wr*[-1 1], n, m);
    L = sparse(1:n, 1:n, w) * (Dr - Dl);
end