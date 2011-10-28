function sol = solve(L, f, u, w)

n = numel(u); % number of basis functions
for k = 1:n
    % v = Lu
    v{k} = L(u{k}); % apply linear operator to each one
end
% L(a_i u_i) = a_i v_i = f
