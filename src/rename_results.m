function rename_results(p)
D = dir(fullfile(p, '*.mat'));

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    src = fullfile(p, f);
    S = load(src);
    Nr = numel(S.g.r)-1;
    Nt = numel(S.g.t)-1;
    Rmax = max(S.g.r);
    t = sprintf('%dx%d_Rmax=%.0f_beta=%.3f_alpha=%.1f_gamma=%.5f.mat', ...
        Nr, Nt, Rmax, S.sol.beta, S.sol.alpha, S.sol.gamma);
    dst = fullfile(p, t);
    cmd = sprintf('mv -v %s %s', src, dst);
    % disp(cmd)
    system(cmd);
    % save(, '-struct', 'S');
end
