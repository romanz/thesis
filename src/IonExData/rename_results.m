D = dir(fullfile('*.mat'));

if isempty(D)
    fprintf('No MAT files found.\n')
    return
end

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    src = fullfile(f);
    S = load(src);
    b = S.sol.beta;
    t = sprintf('%.5f.mat', b);
    dst = fullfile(t);
    cmd = sprintf('mv -v %s %s', src, dst);
    % disp(cmd)
    system(cmd);
    % save(, '-struct', 'S');
end
