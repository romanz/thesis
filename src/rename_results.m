clc; 
clear;
p = '/home/romanz/Research/thesis/src/2012-09-02/e/';
D = dir([p '2012*.mat']);

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    src = [p f];
    S = load(src);
    Nr = numel(S.g.r)-1;
    Nt = numel(S.g.t)-1;
    Rmax = max(S.g.r);
    t = sprintf('%dx%d_Rmax=%.0f_beta=%.3f_alpha=%.1f_Du=%.1f_zeta=%.1f.mat', ...
        Nr, Nt, Rmax, S.sol.beta, S.sol.alpha, S.sol.Du, S.sol.zeta);
    dst = [p t];
    cmd = sprintf('mv -v %s %s', src, dst);
    % disp(cmd)
    system(cmd);
    % save(, '-struct', 'S');
end
