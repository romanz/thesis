clc; 
clear;
p = '1/';
D = dir([p '*.mat']);

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    S = load([p, f]);
    disp(f)
    t = sprintf('%dx%d_Rmax=%.0f_beta=%.3f.mat', ...
        S.Nr-1, S.Nt-1, S.Rmax, S.sol.beta);
    save([p t], '-struct', 'S');
end
