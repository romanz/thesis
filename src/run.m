function run(path, pattern, func)
D = dir(fullfile(path, pattern));
for k = 1:numel(D)
    n = D(k).name;
    func(fullfile(path, n));
end