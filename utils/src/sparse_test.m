function sparse_test()

%     clear sparse_update;
%     mex sparse_update.c;

    n = 1000;
    K = (1:n)';
    I = [K K K];
    J = [K, K+1, K+2];
    sz = [n, n+2];
    tic;
    update = sparse_matrix(I, J, sz);    
    for j=1:10000
        S = update([K -2*K 3*K]);
    end
    toc;
    tic;
    for j=1:10000
        S = sparse(I, J, [K -2*K 3*K], sz(1), sz(2));
    end
    toc;
    return;

    m = 100e3;
    [I, J] = find(sprandn(m, m, 1/m));
    n = numel(I);
    tic;
    for k=1:1000
        updater = sparse_matrix(I, J, [m m]);
        S = updater((1:n)*k);
    end
    toc;
    return

end

function updater = sparse_matrix(I, J, sz)
    assert( numel(I) == numel(J) );
    K = [J(:), I(:), (1:numel(I))'];
    K = sortrows(K);
    K = K(:, 3);
    mask = sparse(I, J, ones(size(I)), sz(1), sz(2));
    function S = update(V)
        assert(numel(V) == numel(I))
        S = mask;
        sparse_update(S, V(K));
    end
    updater = @update;
end
