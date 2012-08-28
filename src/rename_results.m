clc; 
clear;
p = 'better_R_res/';
D = dir([p '*.mat']);

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    src = [p f];
    S = load(src);
    t = sprintf('%dx%d_Rmax=%.0f_beta=%.3f_alpha=%.1f_Du=%.1f_zeta=%.1f.mat', ...
        S.Nr-1, S.Nt-1, S.Rmax, S.sol.beta, S.sol.alpha, S.sol.Du, S.sol.zeta);
    dst = [p t];
    cmd = sprintf('mv -v %s %s', src, dst);
    % disp(cmd)
    system(cmd);
    % save(, '-struct', 'S');
end
