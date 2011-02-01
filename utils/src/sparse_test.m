function sparse_test()

%     clear sparse_update;
%     mex sparse_update.c;

    n = 100000;
    I = 1:n;
    J = 1:n;
    updater = sparse_matrix(I, J, [n n]);
    tic;
    for k=1:1000
        S = updater((1:n)*k);
    end
    toc;
    return

end

function updater = sparse_matrix(I, J, sz)
    I = I(:);
    J = J(:);    
    assert( numel(I) == numel(J) );
    mask = sparse(I, J, ones(size(I)), sz(1), sz(2));
    function S = update(V)
        S = mask;
        sparse_update(S, V);
    end
    updater = @update;
end