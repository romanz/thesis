srcs = {{'sparse_update_inplace.c'}};
for i = 1:numel(srcs)
    clear(srcs{i}{1});
    fprintf('Build %s...', srcs{i}{1})
    mex( srcs{i}{:} );
    fprintf(' OK\n');
end